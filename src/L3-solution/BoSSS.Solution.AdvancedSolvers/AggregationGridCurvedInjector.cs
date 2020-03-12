using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading;

namespace BoSSS.Solution.AdvancedSolvers
{
    public class AggregationGridCurvedInjector
    {
        /// <summary>
        /// The function computes a Np*Np operator for each cell of the grid.
        /// This operator transforms the basis from aggregation level 0 to the
        /// maximum aggregation level. It preserves continuity along inner aggregation cell edges
        /// and the orthonormality of the basis on the aggregation cell. 
        /// Furthermore it seeks a solution that minimizes the gradient jump along inner edges.
        /// </summary>
        /// <param name="_agGrd"></param>
        /// Grid Data for the coarse aggregation level
        /// <param name="_maxDgBasis"></param>
        /// underlying DG basis
        /// <returns> 
        /// - index corresponds to the index of the geometric cell
        /// </returns>
        public static MultidimensionalArray[] AggregateCurvedCells(AggregationGridData _agGrd, Basis _maxDgBasis)
        {
            AggregationGridData agGrd = _agGrd;
            Basis dgBasis = _maxDgBasis;

            // check type of underlying reference elements
            foreach (var refel in agGrd.AncestorGrid.iGeomCells.RefElements)
            {
                if (refel.GetType() != typeof(BoSSS.Foundation.Grid.RefElements.Square)) { throw new NotSupportedException("currently only square elements are supported"); }
            }

            // get dimensions of Injection operator
            int Jagg = agGrd.iLogicalCells.NoOfLocalUpdatedCells;
            int Np = dgBasis.Length;

            MultidimensionalArray[] InjectorCoarseLogicalCell = new MultidimensionalArray[Jagg];
            MultidimensionalArray[] InjectorCoarse = new MultidimensionalArray[agGrd.iGeomCells.Count];

            Console.WriteLine($"Building direct Injector to coarsest level for {Jagg} cells: ");
            var progress = new ProgressBar.ProgressBar();

            // get a list of all inner edges
            int[] iedges;
            List<int> l_iedges = new List<int>();

            for (int i = 0; i < Jagg; i++)
            {
                

                // get parts of current aggregation cell
                int[] parts = agGrd.iLogicalCells.AggregateCellToParts[i];

                l_iedges.Clear();

                for (int j = 0; j < agGrd.iGeomEdges.Count; j++)
                {
                    // check for outer edges and if both cells on the edge belong to the current aggregation cell     
                    if (agGrd.iGeomEdges.LogicalCellIndices[j, 0] == agGrd.iGeomEdges.LogicalCellIndices[j, 1] && agGrd.iGeomEdges.LogicalCellIndices[j, 0] == i)
                    {
                        l_iedges.Add(j);
                    }
                }
                // convert to array
                iedges = l_iedges.ToArray();

                InjectorCoarseLogicalCell[i] = MultidimensionalArray.Create(parts.Length, Np, Np);

                // special case for the first basis function
                InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(-1, 0, 0).Acc(1, GetCoefficients0Degree(agGrd, dgBasis, parts));

                // compute the columns in the Injection operator for subsequent basis functions
                for (int n = 1; n < Np; n++)
                {
                    InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(new int[] { 0, 0, n }, new int[] { parts.Length - 1, n, -1 }).Acc(1, GetCoefficients(agGrd, dgBasis, InjectorCoarseLogicalCell[i], parts, iedges, n));
                }

                // accumulate the direct operator
                foreach(int k in parts)
                {
                    InjectorCoarse[k] = InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(Array.IndexOf(parts, k), -1, -1);
                }

                progress.Report((double)i / Jagg);

            }

            progress.Dispose();
            Console.WriteLine("Done. ");

            return InjectorCoarse;
        }

        // computes the injector coefficients for the first basis functions
        private static MultidimensionalArray GetCoefficients0Degree(AggregationGridData _agGrd, Basis _maxDgBasis, int[] parts)
        {
            var Injector0Degree = MultidimensionalArray.Create(parts.Length);

            // calculate total aggregation cell volume
            double aggVol = 0;
            foreach (int i in parts)
            {
                aggVol += _agGrd.AncestorGrid.Cells.GetCellVolume(i);
            }

            // For the first basis functions with constant value this breaks down to a simple scaling
            double scale = 1 / Math.Sqrt(aggVol);

            int count = 0;
            foreach (int i in parts)
            {
                Injector0Degree[count] = Math.Sqrt(_agGrd.AncestorGrid.Cells.GetCellVolume(i)) * scale;
                count++;
            }

            return Injector0Degree;
        }

        // computes the injector coefficients higher order basis functions
        private static MultidimensionalArray GetCoefficients(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj, int[] parts, int[] iedges, int basisIndex)
        {
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
            aggE1[basisIndex, basisIndex] = 1;
            aggRhs1[basisIndex] = 1;


            // orthogonality
            for (int row = 0; row < basisIndex; row++)
            {
                for (int i = 0; i < partCount; i++)
                {
                    for (int column = 0; column <= row; column++)
                    {
                        aggE1[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
                        //E[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
                    }
                }
            }

            // continuity
            // TODO harder part: Select currentDegree + 1 points per shared edge and ensure continuity there
            // 1. per edge construct currentDegree + 1 points
            // 2. evaluate the basisIndex -basisfunction on the adjacent cells for these points
            // 3. write the values in aggEq
            int rowOffset = basisIndex + 1;
            List<int> usableRowsE = Enumerable.Range(0, basisIndex + 1).ToList();


            for (int row = 0; row < iedges.Length; row++)
            {

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
                if (edgeVertCount < currentDegree + 1)
                {
                    for (int i = 0; i < edgeVertCount; i++)
                    {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }
                    
                    for(int i = edgeVertCount; i < currentDegree + 1; i++)
                    {
                        // construct a point as middle point between the edge vertices
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount, -1));
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount + 1, -1));
                    }
                }
                else
                {
                    for (int i = 0; i < currentDegree + 1; i++)
                    {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }
                }


                NodeSet edge_coords = new NodeSet(refEl, edgeEvalPoints);
                edge_coords.LockForever();

                // evaluate the polynomials on the edge for both cells
                Tuple<MultidimensionalArray, MultidimensionalArray> edge_values = _maxDgBasis.EdgeEval(edge_coords, edge, 1);

                for (int i = 0; i < (currentDegree + 1); i++)
                {

                    for (int column = 0; column < basisIndex + 1; column++)
                    {
                        aggE1[row * (currentDegree + 1) + i + rowOffset, index_i * (basisIndex + 1) + column] = edge_values.Item1[0, i, column];
                        aggE1[row * (currentDegree + 1) + i + rowOffset, index_j * (basisIndex + 1) + column] = -edge_values.Item2[0, i, column];
                        //E[row * (currentDegree + 1) + i + rowOffset, index_i * (basisIndex + 1) + column] = edge_values.Item1[0, i, column];
                        //E[row * (currentDegree + 1) + i + rowOffset, index_j * (basisIndex + 1) + column] = -edge_values.Item2[0, i, column];
                    }
                }

                // test if rows are linear dependent
                usableRowsE.AddRange(CheckLinDependency(aggE1, currentDegree + 1, row * (currentDegree + 1) + rowOffset));                
            }

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
            for (int row = 0; row < iedges.Length; row++)
            {
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
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult)
                    {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        MultidimensionalArray edge_normals = Context.iGeomEdges.NormalsCache.GetNormals_Edge(nodes_ref_edge, edge, 1);
                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_gradient = _maxDgBasis.EdgeEvalGradient(nodes_ref_edge, edge, 1);
                        MultidimensionalArray edge_normal_gradient = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_normal_gradient.ExtractSubArrayShallow(0, -1, -1).Multiply(1.0, edge_gradient.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");
                        edge_normal_gradient.ExtractSubArrayShallow(1, -1, -1).Multiply(-1.0, edge_gradient.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_normal_gradient, edge_normal_gradient, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                                                 {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     gradient.Clear();
                                                     gradient.Acc(1.0, res);
                                                 }).Execute();



                // iterate over both cells (leading coeffcient)
                for (int i = 0; i <= 1; i++)
                {

                    // iterate over all conditions for the current cell
                    for (int k = 1; k < basisIndex + 1; k++)
                    {

                        // iterate over the variables
                        for (int column = 1; column < basisIndex + 1; column++)
                        {
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

            for (int i = 0; i < usableRowsE.Count(); i++)
            {
                E.SetRow(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
                ET.SetColumn(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
            }

            MultidimensionalArray hQ = MultidimensionalArray.Create(usableRowsQ.Count(), varCount);

            for (int i = 0; i < usableRowsQ.Count(); i++)
            {
                hQ.SetRow(i, aggE2.ExtractSubArrayShallow(usableRowsQ.ElementAt(i), -1).To1DArray());
            }

            // solution using quadratic program and sparse matrices
            MultidimensionalArray Q = MultidimensionalArray.Create(hQ.NoOfCols, hQ.NoOfCols);
            double Qscale = 1.0; // 1 / Math.Pow(aggE2.InfNorm(), 2.0);
            Q.DGEMM(Qscale, hQ, hQ, 0.0, true);



            Partitioning qRowPart = new Partitioning(Q.NoOfRows + E.NoOfRows);
            Partitioning qColPart = new Partitioning(Q.NoOfCols + E.NoOfRows);
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

            // restructure to Injector operator form
            for (int i = 0; i < parts.Length; i++)
            {
                InjectorNDegree.ExtractSubArrayShallow(i, -1).Acc(1.0, vector_InjectorNDegree.ExtractSubArrayShallow(new int[] { i * (basisIndex + 1) }, new int[] { (i + 1) * (basisIndex + 1) - 1 }));
            }

            return InjectorNDegree;
        }


        // recursively extracts the Injection operator from ilevel-1 to  ilevel from the direct operator to ilevel
        public static void ExtractInjectorCurved(AggregationGridData[] _agGrd, Basis _maxDgBasis, MultidimensionalArray[][] _Injectors, MultidimensionalArray[] _injectorCoarse, int ilevel)
        {
            Console.WriteLine($"Extracting Injector from level {ilevel - 1} to {ilevel}");
            var progress = new ProgressBar.ProgressBar();

            if (ilevel > 1)
            {

                // get dimensions of Injection operator
                int Jagg = _agGrd[ilevel].iLogicalCells.Count;
                int Np = _maxDgBasis.Length;
                var Context = _maxDgBasis.GridDat;


                _Injectors[ilevel] = new MultidimensionalArray[Jagg];

                MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);
                MultidimensionalArray orthoInv = MultidimensionalArray.Create(Np, Np);

                // compute the mass matrix for lCell on ilevel - 1 when using the Injector to ilevel
                MultidimensionalArray MMtemp = MultidimensionalArray.Create(Np, Np);

                //if (ilevel == _agGrd.Length - 1)
                //{
                //    Console.WriteLine("performing initial orthonormalization");
                //    for (int aggCell = 0; aggCell < Jagg; aggCell++)
                //    {
                //        MMtemp.Clear();
                //        ortho.Clear();
                //        orthoInv.Clear();

                //        int[] parts = _agGrd[ilevel].iLogicalCells.AggregateCellToParts[aggCell];

                //        foreach (int gCell in parts)
                //        {
                //            MMtemp.DGEMM(1.0, _injectorCoarse[gCell], _injectorCoarse[gCell], 1.0, true);
                //        }

                //        MMtemp.AccEye(-1.0);
                //        double dist = MMtemp.L2Norm();
                //        MMtemp.AccEye(1.0);

                //        Console.WriteLine($"intitial deviation in cell {aggCell}: {dist}");

                //        // orthonormalize
                //        // use left or right?
                //        //MMtemp.SymmetricLDLInversion(orthoInv, null);
                //        MMtemp.Cholesky();
                //        MMtemp.TransposeTo(ortho);
                //        ortho.InvertTo(orthoInv);

                //        // form the new direct injector to ilevel - 1
                //        foreach (int gCell in parts)
                //        {
                //            MultidimensionalArray coarseClone = _injectorCoarse[gCell].CloneAs();
                //            _injectorCoarse[gCell].Multiply(1.0, coarseClone, orthoInv, 0.0, "ij", "ik", "kj");
                //        }
                //    }
                //}

                // iterate over all aggregate (logical) cells of ilevel
                for (int aggCell = 0; aggCell < Jagg; aggCell++)
                {

                    // logical cells on ilevel-1 that build the aggregated cell on ilevel
                    int[] aggParts = _agGrd[ilevel].jCellCoarse2jCellFine[aggCell];
                    _Injectors[ilevel][aggCell] = MultidimensionalArray.Create(aggParts.Length, Np, Np);

                    // iterate over all aggregate (logical) cells of ilevel - 1
                    foreach (int lCell in aggParts)
                    {
                        int[] parts = _agGrd[ilevel - 1].iLogicalCells.AggregateCellToParts[lCell];

                        ortho.Clear();
                        orthoInv.Clear();
                        MMtemp.Clear();

                        // based on the cell ordering we have to distinguish different cases
                        // case cell on level k consists of more than one k-1 cells
                        if (aggParts.Length > 1)
                        {

                            // case cell on level k-1 consists of more than one level 0 cells
                            if (parts.Length > 1)
                            {

                                foreach (int gCell in parts)
                                {

                                    //    var cellMask = new CellMask(Context, new[] { new Chunk() { i0 = gCell, Len = 1 } }, MaskType.Geometrical);

                                    //    CellQuadrature.GetQuadrature(new int[2] { Np, Np }, Context,
                                    //    (new CellQuadratureScheme(true, cellMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                                    //    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult)
                                    //    {
                                    //    // get a set of quad points on the reference element
                                    //    NodeSet nodes = QR.Nodes;
                                    //        Debug.Assert(Length == 1);

                                    //        var phi_0 = _maxDgBasis.CellEval(nodes, gCell, 1).ExtractSubArrayShallow(0, -1, -1);
                                    //    // apply the coarse injector
                                    //    MultidimensionalArray phi = phi_0.CloneAs();
                                    //        phi.Multiply(1.0, phi_0, _injectorCoarse[gCell], 0.0, "ki", "kj", "ji");

                                    //        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1);

                                    //        EvalResult.Multiply(1.0, phi, phi, 0.0, "kmn", "kn", "km");
                                    //    },
                                    //    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                                    //                                 {
                                    //                                     Debug.Assert(Length == 1);

                                    //                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1);
                                    //                                     MMtemp.Acc(1.0, res);
                                    //                                 }).Execute();

                                    MMtemp.DGEMM(1.0, _injectorCoarse[gCell], _injectorCoarse[gCell], 1.0, true);

                                }



                                // orthonormalize
                                // use left or right?
                                //MMtemp.SymmetricLDLInversion(orthoInv, null);
                                MMtemp.Cholesky();
                                MMtemp.TransposeTo(ortho);
                                ortho.InvertTo(orthoInv);

                                // form the new direct injector to ilevel - 1
                                foreach (int gCell in parts)
                                {
                                    MultidimensionalArray coarseClone = _injectorCoarse[gCell].CloneAs();
                                    _injectorCoarse[gCell].Multiply(1.0, coarseClone, orthoInv, 0.0, "ij", "ik", "kj");
                                }
                            }
                            // case cell on level k-1 consists of one level 0 cell
                            else
                            {
                                ortho.Acc(1.0, _injectorCoarse[parts[0]]);
                                _injectorCoarse[parts[0]].Clear();
                                _injectorCoarse[parts[0]].AccEye(1.0);
                            }
                        }
                        // case cell on level k consists of one k-1 cell
                        else
                        {
                            ortho.AccEye(1.0);
                        }
                        
                        //orthoInv.InvertTo(ortho);
                        _Injectors[ilevel][aggCell].ExtractSubArrayShallow(Array.IndexOf(aggParts, lCell), -1, -1).Acc(1.0, ortho);



#if DEBUG
                        MMtemp.Clear();

                        // check orthonormality (move to a more suitable location later on)
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

                        // calculate distance to expected unity matrix
                        MultidimensionalArray Eye = MultidimensionalArray.Create(Np, Np);
                        Eye.AccEye(1.0);
                        Eye.Acc(-1.0, MMtemp);
                        double dist = Eye.L2Norm();
                        //Debug.Assert(dist < 1E-10);
                        Console.ForegroundColor = dist < 1E-3 ? ConsoleColor.White : ConsoleColor.Red;
                        //Console.WriteLine($"distance between Eye and MM for level {ilevel - 1} and lCell {lCell}: {dist}", Console.ForegroundColor);
                        Console.ForegroundColor = ConsoleColor.White;
#endif

                    }
                    progress.Report((double)aggCell/Jagg);
                }

                progress.Dispose();
                Console.WriteLine("Done. ");

                // recursive repeat until ilevel = 1
                ExtractInjectorCurved(_agGrd, _maxDgBasis, _Injectors, _injectorCoarse, ilevel - 1);
            }
            else if (ilevel == 1)
            {
                ExtractInjectorCurvedLv1(_agGrd, _maxDgBasis, _Injectors, _injectorCoarse);

                progress.Dispose();
                Console.WriteLine("Done. ");
            }            
        }

        // recursively extracts the Injection operator from ilevel-1 to  ilevel from the direct operator to ilevel
        public static void ExtractInjectorCurvedLv1(AggregationGridData[] _agGrd, Basis _maxDgBasis, MultidimensionalArray[][] _Injectors, MultidimensionalArray[] _injectorCoarse)
        {
            // get dimensions of Injection operator
            int Jagg = _agGrd[1].iLogicalCells.Count;
            int Np = _maxDgBasis.Length;
            var Context = _maxDgBasis.GridDat;


            _Injectors[1] = new MultidimensionalArray[Jagg];


            // iterate over all aggregate (logical) cells of ilevel
            for (int aggCell = 0; aggCell < Jagg; aggCell++)
            {

                // logical cells on ilevel-1 that build the aggregated cell on ilevel
                int[] aggParts = _agGrd[1].jCellCoarse2jCellFine[aggCell];
                _Injectors[1][aggCell] = MultidimensionalArray.Create(aggParts.Length, Np, Np);

                // iterate over all aggregate (logical) cells of ilevel - 1
                foreach (int lCell in aggParts)
                {                   
                    // insert the submatrices from the direct injector 0 to 1
                    _Injectors[1][aggCell].ExtractSubArrayShallow(Array.IndexOf(aggParts, lCell), -1, -1).Acc(1.0, _injectorCoarse[lCell]);

                }
            }
        }

        /// returns the inices of rows in Matrix <paramref name="E"/>
        /// that are linearly independent
        private static List<int> CheckLinDependency(MultidimensionalArray E, int rows, int offset)
        {
            int cols = E.NoOfCols;
            List<int> usableRows = new List<int>();

            int[] pivot = IMatrixExtensions.ReducedRowEchelonForm(E.ExtractSubArrayShallow(new int[] { offset, 0 }, new int[] { offset + rows - 1, cols - 1 }).Transpose()).Item2;

            for (int i = 0; i < pivot.Length; i++) pivot[i] += offset;

            usableRows.AddRange(pivot);
            return usableRows;
        }



    }
}

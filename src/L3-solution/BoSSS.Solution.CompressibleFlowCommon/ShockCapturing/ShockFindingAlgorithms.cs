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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using System;
using System.Collections;
using System.Diagnostics;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {
    /// <summary>
    /// Static methods for locating the inflection points of a DG Field
    /// in two dimensions. Can be used to find the shock location. 
    /// </summary>
    public static class ShockFindingAlgorithms {
        public static NodeSet GetLocalNodeSet(GridData gridData, double[] globalPoint, int globalCellIndex) {
            int D = 2;
            MultidimensionalArray GlobalVerticesIn = MultidimensionalArray.CreateWrapper(globalPoint, 1, D);
            MultidimensionalArray LocalVerticesOut = MultidimensionalArray.CreateWrapper(new double[D], 1, D);

            gridData.TransformGlobal2Local(GlobalVerticesIn, LocalVerticesOut, globalCellIndex, NewtonConvergence: null);

            return new NodeSet(gridData.Cells.GetRefElement(globalCellIndex), LocalVerticesOut);
        }

        public static double[] NormalizedGradient(SinglePhaseField field, int localCellIndex, NodeSet nodeSet) {
            // Evaluate gradient
            int length = 1;
            int noOfNodes = 1;
            int D = 2;
            MultidimensionalArray gradient = MultidimensionalArray.Create(length, noOfNodes, D);
            field.EvaluateGradient(localCellIndex, 1, nodeSet, gradient, ResultCellindexOffset: 0, ResultPreScale: 1.0);

            // Normalize gradient
            double[] tmp = gradient.ExtractSubArrayShallow(0, 0, -1).To1DArray();
            tmp.Normalize();
            return tmp;
        }

        public static double SecondDerivative(SinglePhaseField field, int localCellIndex, NodeSet nodeSet) {
            // Evaluate gradient
            int length = 1;
            int noOfNodes = 1;
            int D = 2;
            MultidimensionalArray gradient = MultidimensionalArray.Create(length, noOfNodes, D);
            field.EvaluateGradient(localCellIndex, 1, nodeSet, gradient, ResultCellindexOffset: 0, ResultPreScale: 1.0);

            // Evalaute Hessian matrix
            MultidimensionalArray hessian = MultidimensionalArray.Create(length, noOfNodes, D, D);
            field.EvaluateHessian(localCellIndex, 1, nodeSet, hessian);

            // Compute second derivative along curve
            double g_alpha_alpha = 2 * ((hessian[0, 0, 0, 0] * gradient[0, 0, 0]
                                       + hessian[0, 0, 0, 1] * gradient[0, 0, 1])
                                       * gradient[0, 0, 0]
                                       + (hessian[0, 0, 1, 0] * gradient[0, 0, 0]
                                       + hessian[0, 0, 1, 1] * gradient[0, 0, 1])
                                       * gradient[0, 0, 1]);

            if (g_alpha_alpha == 0.0 || g_alpha_alpha.IsNaN()) {
                throw new NotSupportedException("Second derivative is zero");
            }

            return g_alpha_alpha;
        }

        /// <summary>
        /// Algorithm to find an inflection point along a curve in the direction of the gradient
        /// </summary>
        /// <param name="gridData">The corresponding grid</param>
        /// <param name="field">The DG field, which shall be evalauted</param>
        /// <param name="maxIterations">The maximum number of iterations (user defined)</param>
        /// <param name="threshold">A threshold, e.g. 1e-6 (user defined)</param>
        /// <param name="points">The points along the curve (first entry has to be user defined)</param>
        /// <param name="secondDerivative">The second derivative at each point</param>
        /// <param name="stepSize">The step size between two points along the curve (first entry has to be user defined)</param>
        /// <returns></returns>
        public static int WalkOnCurve(GridData gridData, SinglePhaseField field, int maxIterations, double threshold, MultidimensionalArray points, double[] secondDerivative, double[] stepSize) {
            // Init
            // Current (global) point
            double[] currentPoint = points.ExtractSubArrayShallow(0, -1).To1DArray();

            // Compute global cell index of current point
            gridData.LocatePoint(currentPoint, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess);

            // Compute local node set
            NodeSet nodeSet = GetLocalNodeSet(gridData, currentPoint, (int)GlobalIndex);

            // Get local cell index of current point
            int j0Grd = gridData.CellPartitioning.i0;
            int jLocal = (int)(GlobalIndex - j0Grd);

            // Evaluate the second derivative
            secondDerivative[0] = SecondDerivative(field, jLocal, nodeSet);

            // Set initial step size to 0.5 * h_minGlobal
            // Set the search direction depending on the sign of the curvature
            stepSize[0] = 0.5 * gridData.Cells.h_minGlobal;
            //if (secondDerivative[0] < 0.0) {
            //    stepSize[0] = -stepSize[0];
            //}

            int n = 1;
            while (n < maxIterations + 1) {
                // Evaluate the gradient of the current point
                double[] gradient = NormalizedGradient(field, jLocal, nodeSet);

                // Compute new point along curve
                points[n, 0] = currentPoint[0] + gradient[0] * stepSize[n - 1];
                points[n, 1] = currentPoint[1] + gradient[1] * stepSize[n - 1];
                currentPoint = points.ExtractSubArrayShallow(n, -1).To1DArray();

                // Check if new point is still in the same cell or has moved to one of its neighbours
                if (!gridData.Cells.IsInCell(currentPoint, jLocal)) {
                    // Get CellMask of cell neighbours
                    gridData.GetCellNeighbours(jLocal, GetCellNeighbours_Mode.ViaVertices, out int[] cellNeighbours, out int[] connectingEntities);
                    BitArray btArray = new BitArray(gridData.Cells.NoOfLocalUpdatedCells);
                    foreach (int neighbour in cellNeighbours) {
                        btArray[neighbour] = true;
                    }
                    CellMask neighbours = new CellMask(gridData, btArray);

                    // Compute global cell index of new point
                    gridData.LocatePoint(currentPoint, out GlobalId, out GlobalIndex, out IsInside, out OnThisProcess, neighbours);
                    if (!IsInside) {
                        Console.WriteLine("New point is outside of the grid. Skip this curve and go on.");
                        break;
                    }

                    // Compute local node set
                    nodeSet = GetLocalNodeSet(gridData, currentPoint, (int)GlobalIndex);

                    // Get local cell index of new point
                    jLocal = (int)(GlobalIndex - j0Grd);
                }

                // Evaluate the second derivative of the new point
                NotSupportedException e = null;
                try {
                    secondDerivative[n] = SecondDerivative(field, jLocal, nodeSet);
                } catch (NotSupportedException ee) {
                    e = ee;
                    break;
                }

                // Check if sign of second derivative has changed
                // Halve the step size and change direction
                if (Math.Sign(secondDerivative[n]) != Math.Sign(secondDerivative[n - 1])) {
                    stepSize[n] = -stepSize[n - 1] / 2;
                } else {
                    stepSize[n] = stepSize[n - 1];
                }

                n++;

                // Termination criterion
                // Remark: n has already been incremented!
                double[] diff = new double[currentPoint.Length];
                diff[0] = points[n - 1, 0] - points[n - 2, 0];
                diff[1] = points[n - 1, 1] - points[n - 2, 1];
                if (diff.L2Norm() < threshold) {
                    break;
                }
            }

            return n;
        }

        /// <summary>
        /// Returns all cells with the respective artificial viscosity value
        /// if the value is larger than zero
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="avField"></param>
        /// <returns><see cref="MultidimensionalArray"/>, first index: cellIndex, second index: artificial viscosity value</returns>
        public static MultidimensionalArray GetAvMeanValues(GridData gridData, SinglePhaseField avField) {
            CellMask allCells = CellMask.GetFullMask(gridData);
            double[] cellIndices = new double[allCells.NoOfItemsLocally];
            double[] avValues = new double[allCells.NoOfItemsLocally];

            int count = 0;
            foreach (int cell in allCells.ItemEnum) {
                double avValue = avField.GetMeanValue(cell);
                if (avValue > 0.0) {
                    cellIndices[count] = cell;
                    avValues[count] = avValue;
                    count++;
                }
            }

            Array.Resize(ref cellIndices, count);
            Array.Resize(ref avValues, count);

            MultidimensionalArray result = MultidimensionalArray.Create(count, 2);
            result.SetSubVector(cellIndices, new int[] { -1, 0 });
            result.SetSubVector(avValues, new int[] { -1, 1 });

            return result;
        }

        /// <summary>
        /// Get a 3x3 grid in a cell
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="jCell">Local cell index</param>
        /// <returns><see cref="MultidimensionalArray"/>, first index: x-coord., second index: y-coord.</returns>
        public static MultidimensionalArray GetGlobal3x3NodeSet(GridData gridData, int jCell) {
            int D = 2;
            int numOfNodes = 9;
            double hmin = gridData.Cells.h_minGlobal;
            double[] cellCenter = gridData.Cells.GetCenter(jCell);

            MultidimensionalArray grid = MultidimensionalArray.Create(numOfNodes, D);

            // Create a 3x3 grid in a cell
            // Bottom row
            // Node 0
            grid[0, 0] = cellCenter[0] - 0.4 * hmin;
            grid[0, 1] = cellCenter[1] - 0.4 * hmin;

            // Node 1
            grid[1, 0] = cellCenter[0];
            grid[1, 1] = cellCenter[1] - 0.4 * hmin;

            // Node 2
            grid[2, 0] = cellCenter[0] + 0.4 * hmin;
            grid[2, 1] = cellCenter[1] - 0.4 * hmin;

            // Middle row
            // Node 3
            grid[3, 0] = cellCenter[0] - 0.4 * hmin;
            grid[3, 1] = cellCenter[1];

            // Node 4
            grid[4, 0] = cellCenter[0];
            grid[4, 1] = cellCenter[1];

            // Node 5
            grid[5, 0] = cellCenter[0] + 0.4 * hmin;
            grid[5, 1] = cellCenter[1];

            // Top row
            // Node 6
            grid[6, 0] = cellCenter[0] - 0.4 * hmin;
            grid[6, 1] = cellCenter[1] + 0.4 * hmin;

            // Node 7
            grid[7, 0] = cellCenter[0];
            grid[7, 1] = cellCenter[1] + 0.4 * hmin;

            // Node 8
            grid[8, 0] = cellCenter[0] + 0.4 * hmin;
            grid[8, 1] = cellCenter[1] + 0.4 * hmin;

            return grid;
        }

        public static MultidimensionalArray SortPoints(GridData gridData, SinglePhaseField field, MultidimensionalArray points) {
            double[] sortedXCoords = new double[points.Lengths[0]];
            double[] sortedYCoords = new double[points.Lengths[0]];
            int keptEntries = 0;

            for (int i = 0; i < points.Lengths[0]; i++) {
                // Compute global cell index of the point
                double[] currentPoint = points.ExtractSubArrayShallow(i, 0, -1).To1DArray();
                gridData.LocatePoint(currentPoint, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess);

                // Compute local node set
                NodeSet nodeSet = GetLocalNodeSet(gridData, currentPoint, (int)GlobalIndex);

                // Get local cell index of current point
                int j0Grd = gridData.CellPartitioning.i0;
                int jLocal = (int)(GlobalIndex - j0Grd);

                // Calculate secondDerivative
                double secondDerivative = SecondDerivative(field, jLocal, nodeSet);

                // Select points where second derivative is larger than zero
                if (secondDerivative > 0.0) {
                    sortedXCoords[keptEntries] = currentPoint[0];
                    sortedYCoords[keptEntries] = currentPoint[1];
                    keptEntries++;
                }
            }

            // Resize final array
            MultidimensionalArray result = MultidimensionalArray.Create(keptEntries, points.Lengths[1], points.Lengths[2]);
            for (int i = 0; i < keptEntries; i++) {
                result[i, 0, 0] = sortedXCoords[i];
                result[i, 0, 1] = sortedYCoords[i];
            }

            return result;
        }
    }
}

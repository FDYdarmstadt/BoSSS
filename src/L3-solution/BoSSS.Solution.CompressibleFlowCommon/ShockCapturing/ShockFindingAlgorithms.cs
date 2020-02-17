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
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
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

        private static SinglePhaseField gradientX;
        private static SinglePhaseField gradientY;
        private static SinglePhaseField hessianXX;
        private static SinglePhaseField hessianXY;
        private static SinglePhaseField hessianYX;
        private static SinglePhaseField hessianYY;

        private static MultidimensionalArray resultGradX;
        private static MultidimensionalArray resultGradY;
        private static MultidimensionalArray resultHessXX;
        private static MultidimensionalArray resultHessXY;
        private static MultidimensionalArray resultHessYX;
        private static MultidimensionalArray resultHessYY;

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

        public static double[] NormalizedGradientByFlux(SinglePhaseField field, int jCell, NodeSet nodeSet) {
            if (gradientX == null) {
                // Evaluate gradient
                Basis basis = new Basis(field.GridDat, field.Basis.Degree);
                gradientX = new SinglePhaseField(basis, "gradientX");
                gradientY = new SinglePhaseField(basis, "gradientY");
                gradientX.DerivativeByFlux(1.0, field, d: 0);
                gradientY.DerivativeByFlux(1.0, field, d: 1);
            }

            if (resultGradX == null) {
                resultGradX = MultidimensionalArray.Create(1, 1);
                resultGradY = MultidimensionalArray.Create(1, 1);
            }

            gradientX.Evaluate(jCell, 1, nodeSet, resultGradX);
            gradientY.Evaluate(jCell, 1, nodeSet, resultGradY);

            double[] gradient = new double[] { resultGradX[0, 0], resultGradY[0, 0] };
            gradient.Normalize();
            return gradient;
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

        public static double SecondDerivativeByFlux(SinglePhaseField field, int jCell, NodeSet nodeSet) {
            if (gradientX == null) {
                // Evaluate gradient
                Basis basis = new Basis(field.GridDat, field.Basis.Degree);
                gradientX = new SinglePhaseField(basis, "gradientX");
                gradientY = new SinglePhaseField(basis, "gradientY");
                gradientX.DerivativeByFlux(1.0, field, d: 0);
                gradientY.DerivativeByFlux(1.0, field, d: 1);
            }

            if (hessianXX == null) {
                // Evaluate Hessian matrix
                Basis basis = new Basis(field.GridDat, field.Basis.Degree);
                hessianXX = new SinglePhaseField(basis, "hessianXX");
                hessianXY = new SinglePhaseField(basis, "hessianXY");
                hessianYX = new SinglePhaseField(basis, "hessianYX");
                hessianYY = new SinglePhaseField(basis, "hessianYY");
                hessianXX.DerivativeByFlux(1.0, gradientX, d: 0);
                hessianXY.DerivativeByFlux(1.0, gradientX, d: 1);
                hessianYX.DerivativeByFlux(1.0, gradientY, d: 0);
                hessianYY.DerivativeByFlux(1.0, gradientY, d: 1);

                if (PatchRecoveryOn) {
                    hessianXX = PatchRecovery(hessianXX);
                    hessianXY = PatchRecovery(hessianXY);
                    hessianYX = PatchRecovery(hessianYX);
                    hessianYY = PatchRecovery(hessianYY);
                }
            }

            if (resultGradX == null) {
                resultGradX = MultidimensionalArray.Create(1, 1);
                resultGradY = MultidimensionalArray.Create(1, 1);
            }

            if (resultHessXX == null) {
                resultHessXX = MultidimensionalArray.Create(1, 1);
                resultHessXY = MultidimensionalArray.Create(1, 1);
                resultHessYX = MultidimensionalArray.Create(1, 1);
                resultHessYY = MultidimensionalArray.Create(1, 1);
            }

            gradientX.Evaluate(jCell, 1, nodeSet, resultGradX);
            gradientY.Evaluate(jCell, 1, nodeSet, resultGradY);

            // Use patch recovery if it has been activated in method "WalkOnCurve(...)"


            hessianXX.Evaluate(jCell, 1, nodeSet, resultHessXX);
            hessianXY.Evaluate(jCell, 1, nodeSet, resultHessXY);
            hessianYX.Evaluate(jCell, 1, nodeSet, resultHessYX);
            hessianYY.Evaluate(jCell, 1, nodeSet, resultHessYY);

            // Compute second derivative along curve
            double g_alpha_alpha = 2 * ((resultHessXX[0, 0] * resultGradX[0, 0] + resultHessXY[0, 0] * resultGradY[0, 0]) * resultGradX[0, 0]
                + (resultHessYX[0, 0] * resultGradX[0, 0] + resultHessYY[0, 0] * resultGradY[0, 0]) * resultGradY[0, 0]);

            if (g_alpha_alpha == 0.0 || g_alpha_alpha.IsNaN()) {
                throw new NotSupportedException("Second derivative is zero");
            }

            return g_alpha_alpha;
        }

        private static bool PatchRecoveryOn;

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
        public static int WalkOnCurve(GridData gridData, SinglePhaseField field, int maxIterations, double threshold, MultidimensionalArray points, double[] secondDerivative, double[] stepSize, double[] values, out bool converged, bool patchRecoveryOn = false, bool byFlux = true) {
            // Init
            converged = false;
            PatchRecoveryOn = patchRecoveryOn;

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
            //NotSupportedException e = null;
            try {
                if (byFlux) {
                    secondDerivative[0] = SecondDerivativeByFlux(field, jLocal, nodeSet);
                } else {
                    secondDerivative[0] = SecondDerivative(field, jLocal, nodeSet);
                }
            } catch (NotSupportedException ee) {
                //e = ee;
                return 0;
            }

            // Evaluate the function
            MultidimensionalArray f = MultidimensionalArray.Create(1, 1);
            field.Evaluate(jLocal, 1, nodeSet, f);
            values[0] = f[0, 0];

            // Set initial step size to 0.5 * h_minGlobal
            // Set the search direction depending on the sign of the curvature
            stepSize[0] = 0.5 * gridData.Cells.h_minGlobal;

            int n = 1;
            while (n < maxIterations + 1) {
                // Evaluate the gradient of the current point
                double[] gradient;
                if (byFlux) {
                    gradient = NormalizedGradientByFlux(field, jLocal, nodeSet);
                } else {
                    gradient = NormalizedGradient(field, jLocal, nodeSet);
                }

                // Compute new point along curve
                currentPoint[0] = currentPoint[0] + gradient[0] * stepSize[n - 1];
                currentPoint[1] = currentPoint[1] + gradient[1] * stepSize[n - 1];

                // New point has been calculated --> old node set is invalid
                nodeSet = null;

                // Check if new point is still in the same cell or has moved to one of its neighbours
                if (!gridData.Cells.IsInCell(currentPoint, jLocal)) {
                    // Get indices of cell neighbours
                    gridData.GetCellNeighbours(jLocal, GetCellNeighbours_Mode.ViaVertices, out int[] cellNeighbours, out int[] connectingEntities);

                    double[] newLocalCoord = new double[currentPoint.Length];
                    bool found = false;
                    // Find neighbour
                    foreach (int neighbour in cellNeighbours) {
                        if (gridData.Cells.IsInCell(currentPoint, neighbour, newLocalCoord)) {
                            // If neighbour has been found, update
                            jLocal = neighbour + j0Grd;
                            nodeSet = GetLocalNodeSet(gridData, currentPoint, jLocal);
                            found = true;
                            break;
                        }
                    }

                    if (found == false) {
                        Console.WriteLine("No neighbour found");
                        return n;
                    }
                } else {
                    // New point is still in the same cell --> update only the coordiantes (local node set)
                    nodeSet = GetLocalNodeSet(gridData, currentPoint, jLocal + j0Grd);
                }

                // Update output
                points[n, 0] = currentPoint[0];
                points[n, 1] = currentPoint[1];

                // Evaluate the function
                f.Clear();
                field.Evaluate(jLocal, 1, nodeSet, f);
                values[n] = f[0, 0];

                // Evaluate the second derivative of the new point
                //NotSupportedException e = null;
                try {
                    if (byFlux) {
                        secondDerivative[n] = SecondDerivativeByFlux(field, jLocal, nodeSet);
                    } else {
                        secondDerivative[n] = SecondDerivative(field, jLocal, nodeSet);
                    }
                } catch (NotSupportedException ee) {
                    //e = ee;
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
                    converged = true;
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
        /// Get a grid consisting of 3 x 3 points per cell
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="jCell">Local cell index</param>
        /// <returns><see cref="MultidimensionalArray"/>, first index: x-coord., second index: y-coord.</returns>
        public static MultidimensionalArray Get3x3Grid(GridData gridData, int jCell) {
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

        public static MultidimensionalArray SortPoints(GridData gridData, SinglePhaseField field, MultidimensionalArray points, bool byFlux = false) {
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
                double secondDerivative;
                if (byFlux) {
                    secondDerivative = SecondDerivativeByFlux(field, jLocal, nodeSet);
                } else {
                    secondDerivative = SecondDerivative(field, jLocal, nodeSet);
                }

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

        /// <summary>
        /// Performs the patch recovery on a given DG field
        /// </summary>
        /// <param name="input">A <see cref="SinglePhaseField"/></param>
        /// <returns>A patch recovered <see cref="SinglePhaseField"/> with a degree of inputDegreee + 2</returns>
        public static SinglePhaseField PatchRecovery(SinglePhaseField input) {
            Console.WriteLine(String.Format("Patch recovery of field {0} started", input.Identification));

            Basis prcBasis = new Basis(input.GridDat, input.Basis.Degree + 2);
            L2PatchRecovery prc = new L2PatchRecovery(input.Basis, prcBasis, CellMask.GetFullMask(input.GridDat), RestrictToCellMask: true);
            SinglePhaseField prcField = new SinglePhaseField(prcBasis);
            prc.Perform(prcField, input);

            Console.WriteLine(String.Format("Patch recovery of field {0} finished", input.Identification));

            return prcField;
        }
    }
}

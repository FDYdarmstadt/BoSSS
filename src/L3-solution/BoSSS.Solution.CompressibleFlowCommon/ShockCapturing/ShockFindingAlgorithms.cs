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

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {
    /// <summary>
    /// Static methods for locating the inflection points of a DG Field
    /// in two dimensions. Can be used to find the shock location. 
    /// </summary>
    public static class ShockFindingAlgorithms {
        public static NodeSet GetLocalNodeSet(GridData gridData, double[] globalPoint, int globalIndex) {
            int D = 2;
            MultidimensionalArray GlobalVerticesIn = MultidimensionalArray.CreateWrapper(globalPoint, 1, D);
            MultidimensionalArray LocalVerticesOut = MultidimensionalArray.CreateWrapper(new double[D], 1, D);

            gridData.TransformGlobal2Local(GlobalVerticesIn, LocalVerticesOut, globalIndex, NewtonConvergence: null);

            return new NodeSet(gridData.Cells.GetRefElement(globalIndex), LocalVerticesOut);
        }

        public static double[] Gradient(SinglePhaseField field, int localIndex, NodeSet nodeSet) {
            // Evaluate gradient
            int length = 1;
            int noOfNodes = 1;
            int D = 2;
            MultidimensionalArray gradient = MultidimensionalArray.Create(length, noOfNodes, D);
            field.EvaluateGradient(localIndex, 1, nodeSet, gradient, ResultCellindexOffset: 0, ResultPreScale: 1.0);

            return gradient.ExtractSubArrayShallow(0, 0, -1).To1DArray();
        }

        public static double SecondDerivative(SinglePhaseField field, int localIndex, NodeSet nodeSet) {
            // Evaluate gradient
            int length = 1;
            int noOfNodes = 1;
            int D = 2;
            MultidimensionalArray gradient = MultidimensionalArray.Create(length, noOfNodes, D);
            field.EvaluateGradient(localIndex, 1, nodeSet, gradient, ResultCellindexOffset: 0, ResultPreScale: 1.0);

            // Evalaute Hessian matrix
            MultidimensionalArray hessian = MultidimensionalArray.Create(length, noOfNodes, D, D);
            field.EvaluateHessian(localIndex, 1, nodeSet, hessian);

            // Compute second derivative along curve
            double g_alpha_alpha = 2 * ((hessian[0, 0, 0, 0] * gradient[0, 0, 0]
                                       + hessian[0, 0, 0, 1] * gradient[0, 0, 1])
                                       * gradient[0, 0, 0]
                                       + (hessian[0, 0, 1, 0] * gradient[0, 0, 0]
                                       + hessian[0, 0, 1, 1] * gradient[0, 0, 1])
                                       * gradient[0, 0, 1]);

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
            double[] initialPoint = points.ExtractSubArrayShallow(0, -1).To1DArray();

            // Compute global cell index of current point
            gridData.LocatePoint(initialPoint, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess);

            // Compute local node set
            NodeSet nodeSet = GetLocalNodeSet(gridData, initialPoint, (int)GlobalIndex);

            // Get local cell index of current point
            int j0Grd = gridData.CellPartitioning.i0;
            int jLocal = (int)(GlobalIndex - j0Grd);

            secondDerivative[0] = SecondDerivative(field, jLocal, nodeSet);

            int n = 1;
            while (n < maxIterations + 1) {
                // Current (global) point
                double[] currentPoint = points.ExtractSubArrayShallow(n - 1, -1).To1DArray();

                // Compute global cell index of current point
                gridData.LocatePoint(currentPoint, out GlobalId, out GlobalIndex, out IsInside, out OnThisProcess);

                // Compute local node set
                nodeSet = GetLocalNodeSet(gridData, currentPoint, (int)GlobalIndex);

                // Get local cell index of current point
                j0Grd = gridData.CellPartitioning.i0;
                jLocal = (int)(GlobalIndex - j0Grd);

                // Evaluate the gradient of the current point
                double[] gradient = Gradient(field, jLocal, nodeSet);

                //gridData.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaVertices, out int[] cellneighbours, out int[] connectingEntities);
                //gridData.Cells.IsInCell()       

                // Compute new point along curve
                double[] newPoint = new double[currentPoint.Length];
                newPoint[0] = currentPoint[0] + gradient[0] * stepSize[n - 1];
                newPoint[1] = currentPoint[1] + gradient[1] * stepSize[n - 1];
                points[n, 0] = newPoint[0];
                points[n, 1] = newPoint[1];
                               
                // Compute global cell index of new point
                gridData.LocatePoint(newPoint, out GlobalId, out GlobalIndex, out IsInside, out OnThisProcess);

                // Compute local node set
                nodeSet = GetLocalNodeSet(gridData, newPoint, (int)GlobalIndex);

                // Get local cell index of current point
                jLocal = (int)(GlobalIndex - j0Grd);

                // Evaluate the second derivative of the new point
                secondDerivative[n] = SecondDerivative(field, jLocal, nodeSet);

                // Check if sign of second derivative has changed
                // Halve the step size and change direction
                if (Math.Sign(secondDerivative[n]) != Math.Sign(secondDerivative[n - 1])) {
                    stepSize[n] = -stepSize[n - 1] / 2;
                } else {
                    stepSize[n] = stepSize[n - 1];
                }

                n++;

                // Termination criterion
                double[] diff = new double[newPoint.Length];
                diff[0] = newPoint[0] - currentPoint[0];
                diff[1] = newPoint[1] - currentPoint[1];
                if (diff.L2Norm() < threshold) {
                    break;
                }
            }

            return n;
        }

        public static double[] GetAVMeanValues(GridData gridData, SinglePhaseField avField) {
            CellMask allCells = CellMask.GetFullMask(gridData);
            int numOfCells = gridData.iGeomCells.Count;
            double[] result = new double[numOfCells];

            foreach (int cell in allCells.ItemEnum) {
                result[cell] = avField.GetMeanValue(cell);
            }

            return result;
        }
    }
}

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
using BoSSS.Foundation.IO;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {
    public class InflectionPointFinder {
        private SinglePhaseField gradientX;
        private SinglePhaseField gradientY;
        private SinglePhaseField hessianXX;
        private SinglePhaseField hessianXY;
        private SinglePhaseField hessianYX;
        private SinglePhaseField hessianYY;

        private MultidimensionalArray resultGradX;
        private MultidimensionalArray resultGradY;
        private MultidimensionalArray resultHessXX;
        private MultidimensionalArray resultHessXY;
        private MultidimensionalArray resultHessYX;
        private MultidimensionalArray resultHessYY;

        private readonly GridData gridData;
        private readonly SinglePhaseField densityField;
        private readonly SinglePhaseField avField;
        private readonly SinglePhaseField levelSetField;

        private readonly string sessionPath;
        private readonly ISessionInfo session;

        private bool patchRecoveryGradient;
        private bool patchRecoveryHessian;

        /// <summary>
        /// Main result data structure
        /// <see cref="MultidimensionalArray"/> with three dimensions
        /// [0]: number of seedings points, [1]: maximum number of iterations to find the inflection point, [2]: data (*)
        /// Explanation for (*):
        /// [0]: x-coordinates, [1]: y-coordinates, [2]: field values, [3]: second derivative, [4]: step size
        /// </summary>
        public MultidimensionalArray MorePoints {
            get;
            private set;
        }

        bool[] converged;
        int[] iterations;

        //####################################################################################################
        public InflectionPointFinder(string sessionPath, ISessionInfo session) {
            this.sessionPath = sessionPath;
            this.session = session;

            ITimestepInfo myTimestep = session.Timesteps.Last();
            this.gridData = (GridData)myTimestep.Fields.First().GridDat;
            this.densityField = (SinglePhaseField)myTimestep.Fields.Find("rho");
            this.avField = (SinglePhaseField)myTimestep.Fields.Find("artificialViscosity");
            this.levelSetField = (SinglePhaseField)myTimestep.Fields.Find("levelSet");
        }

        public MultidimensionalArray FindPoints(SeedingSetup testCase, bool patchRecoveryGradient, bool patchRecoveryHessian, int maxNumOfIterations = 100, double eps = 1e-12) {
            #region Init
            this.patchRecoveryGradient = patchRecoveryGradient;
            this.patchRecoveryHessian = patchRecoveryHessian;
            #endregion

            int numOfPoints;
            int[] jCell;

            #region Create seedings points based on artificial viscosity
            MultidimensionalArray avValues = ShockFindingHelperFunctions.GetAVMeanValues(gridData, avField);

            switch (testCase) {
                case SeedingSetup.av:
                    numOfPoints = avValues.Lengths[0];
                    MorePoints = MultidimensionalArray.Create(numOfPoints, maxNumOfIterations + 1, 5);
                    iterations = new int[numOfPoints];
                    converged = new bool[numOfPoints];
                    jCell = new int[numOfPoints];

                    // Seed points
                    for (int i = 0; i < numOfPoints; i++) {
                        int jLocal = (int)avValues[i, 0];
                        double[] cellCenter = gridData.Cells.GetCenter(jLocal);
                        MorePoints[i, 0, 0] = cellCenter[0];
                        MorePoints[i, 0, 1] = cellCenter[1];
                    }
                    break;

                case SeedingSetup.av3x3:
                    int numOfCells = avValues.Lengths[0];
                    int numOfPointsPerCell = 9;
                    numOfPoints = numOfCells * numOfPointsPerCell;
                    MorePoints = MultidimensionalArray.Create(numOfPoints, maxNumOfIterations + 1, 5);
                    iterations = new int[numOfPoints];
                    converged = new bool[numOfPoints];
                    jCell = new int[numOfPoints];

                    // Seed points
                    for (int i = 0; i < numOfCells; i++) {
                        int jLocal = (int)avValues[i, 0];
                        MultidimensionalArray grid = ShockFindingHelperFunctions.Get3x3Grid(gridData, jLocal);
                        for (int j = 0; j < numOfPointsPerCell; j++) {
                            MorePoints[i * numOfPointsPerCell + j, 0, 0] = grid[j, 0];
                            MorePoints[i * numOfPointsPerCell + j, 0, 1] = grid[j, 1];
                        }
                    }
                    break;

                default:
                    throw new NotSupportedException("This setting does not exist.");
            }
            Console.WriteLine("Total number of seeding points: " + MorePoints.Lengths[0]);
            #endregion

            #region Find inflection point for every seeding point
            double[] x;
            double[] y;

            for (int i = 0; i < MorePoints.Lengths[0]; i++) {
                /*
                MultidimensionalArray test = MultidimensionalArray.Create(10, 2, 5);
                test[0, 0, 0] = 0;
                test[0, 0, 1] = 1;
                test[0, 0, 2] = 2;
                test[0, 0, 3] = 3;
                test[0, 1, 4] = 4;
                test[0, 1, 1] = 1;
                test[0, 1, 2] = 2;
                test[0, 1, 3] = 3;
                test[0, 1, 4] = 4;
                var bla = test.ExtractSubArrayShallow(0, -1, -1);
                var bla3 = test.ExtractSubArrayShallow(0, -1, 3);
                */

                int iter = WalkOnCurve(gridData, densityField, MorePoints.ExtractSubArrayShallow(i, -1, -1), out bool pointFound, out int jLocal);

                x = MorePoints.ExtractSubArrayShallow(i, -1, 0).To1DArray();
                y = MorePoints.ExtractSubArrayShallow(i, -1, 1).To1DArray();
                Array.Resize(ref x, iter);
                Array.Resize(ref y, iter);
                //double[] results = new double[iter];

                // NECESSARY????????
                //switch (testCase) {
                //    case SeedingSetup.av:
                //    case SeedingSetup.av3x3:
                //        //fValues[i, j] = field.ProbeAt(new double[] {x[j], y[j]});
                //        fValues.SetSubVector(values, new int[] { i, -1 });
                //        break;
                //}

                // Save more stuff
                iterations[i] = iter;
                converged[i] = pointFound;
                jCell[i] = jLocal;

                if (i % 100 == 0) {
                    Console.WriteLine("Point " + i);
                }
            }
            #endregion

            return MorePoints;
        }

        public void PlotAction(bool plotDGFields, bool plotSeedingsPoints, bool plotInflectionsPoints, bool plotCurves, bool plotStartEndPairs) {
            if (plotDGFields) {
                Tecplot.Tecplot plotDriver = new Tecplot.Tecplot(gridData, showJumps: true, ghostZone: false, superSampling: 2);
                plotDriver.PlotFields(sessionPath + "Fields", 0.0, new DGField[] { densityField, avField, levelSetField });
            }

            if (plotSeedingsPoints) {
                // Write text file
                using (System.IO.StreamWriter sw = new System.IO.StreamWriter(sessionPath + "seedingPoints.txt")) {
                    string resultLine;
                    for (int j = 0; j < MorePoints.Lengths[0]; j++) {
                        resultLine = MorePoints[j, 0, 0] + "\t"
                                   + MorePoints[j, 0, 1] + "\t"
                                   + MorePoints[j, 0, 2];
                        sw.WriteLine(resultLine);
                    }
                    sw.Flush();
                }
            }

            if (plotInflectionsPoints) {
                // Write text file
                using (System.IO.StreamWriter sw = new System.IO.StreamWriter(sessionPath + "inflectionPoints.txt")) {
                    string resultLine;
                    for (int j = 0; j < MorePoints.Lengths[0]; j++) {
                        int pointFound = converged[j] ? 1 : 0;
                        resultLine = MorePoints[j, iterations[j] - 1, 0] + "\t"
                                   + MorePoints[j, iterations[j] - 1, 1] + "\t"
                                   + MorePoints[j, iterations[j] - 1, 2] + "\t"
                                   + pointFound;
                        sw.WriteLine(resultLine);
                    }

                    sw.Flush();
                }
            }

            if (plotCurves) {
                for (int i = 0; i < MorePoints.Lengths[0]; i++) {
                    using (System.IO.StreamWriter sw = new System.IO.StreamWriter(sessionPath + "curve_" + i + ".txt")) {
                        string resultLine;
                        for (int j = 0; j < iterations[i]; j++) {
                            resultLine = MorePoints[i, j, 0] + "\t" + MorePoints[i, j, 1] + "\t" + MorePoints[i, j, 2];
                            sw.WriteLine(resultLine);
                        }
                        sw.Flush();
                    }
                }
            }

            for (int j = 0; j < MorePoints.Lengths[0]; j++) {
                using (System.IO.StreamWriter sw = new System.IO.StreamWriter(sessionPath + "startEndPairs_" + j + ".txt")) {

                    string resultLine;
                    int pointFound = converged[j] ? 1 : 0;

                    // Starting point
                    resultLine = MorePoints[j, 0, 0] + "\t"
                               + MorePoints[j, 0, 1] + "\t"
                               + MorePoints[j, 0, 2] + "\t"
                               + pointFound;
                    sw.WriteLine(resultLine);

                    // End point
                    resultLine = MorePoints[j, iterations[j] - 1, 0] + "\t"
                               + MorePoints[j, iterations[j] - 1, 1] + "\t"
                               + MorePoints[j, iterations[j] - 1, 2] + "\t"
                               + pointFound;
                    sw.WriteLine(resultLine);

                    sw.Flush();
                }
            }
        }

        //####################################################################################################

        private double[] NormalizedGradientByFlux(SinglePhaseField field, int jCell, NodeSet nodeSet) {
            if (this.gradientX == null) {
                // Evaluate gradient
                Basis basis = new Basis(field.GridDat, field.Basis.Degree);
                gradientX = new SinglePhaseField(basis, "gradientX");
                gradientY = new SinglePhaseField(basis, "gradientY");
                gradientX.DerivativeByFlux(1.0, field, d: 0);
                gradientY.DerivativeByFlux(1.0, field, d: 1);

                if (this.patchRecoveryGradient) {
                    gradientX = ShockFindingHelperFunctions.PatchRecovery(gradientX);
                    gradientY = ShockFindingHelperFunctions.PatchRecovery(gradientY);
                }
            }

            if (this.resultGradX == null) {
                resultGradX = MultidimensionalArray.Create(1, 1);
                resultGradY = MultidimensionalArray.Create(1, 1);
            }

            gradientX.Evaluate(jCell, 1, nodeSet, resultGradX);
            gradientY.Evaluate(jCell, 1, nodeSet, resultGradY);

            double[] gradient = new double[] { resultGradX[0, 0], resultGradY[0, 0] };
            gradient.Normalize();
            return gradient;
        }

        private double SecondDerivativeByFlux(SinglePhaseField field, int jCell, NodeSet nodeSet) {
            if (this.gradientX == null) {
                // Evaluate gradient
                Basis basis = new Basis(field.GridDat, field.Basis.Degree);
                gradientX = new SinglePhaseField(basis, "gradientX");
                gradientY = new SinglePhaseField(basis, "gradientY");
                gradientX.DerivativeByFlux(1.0, field, d: 0);
                gradientY.DerivativeByFlux(1.0, field, d: 1);

                if (this.patchRecoveryGradient) {
                    gradientX = ShockFindingHelperFunctions.PatchRecovery(gradientX);
                    gradientY = ShockFindingHelperFunctions.PatchRecovery(gradientY);
                }
            }

            if (this.hessianXX == null) {
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

                if (this.patchRecoveryHessian) {
                    hessianXX = ShockFindingHelperFunctions.PatchRecovery(hessianXX);
                    hessianXY = ShockFindingHelperFunctions.PatchRecovery(hessianXY);
                    hessianYX = ShockFindingHelperFunctions.PatchRecovery(hessianYX);
                    hessianYY = ShockFindingHelperFunctions.PatchRecovery(hessianYY);
                }
            }

            if (this.resultGradX == null) {
                resultGradX = MultidimensionalArray.Create(1, 1);
                resultGradY = MultidimensionalArray.Create(1, 1);
            }

            if (this.resultHessXX == null) {
                resultHessXX = MultidimensionalArray.Create(1, 1);
                resultHessXY = MultidimensionalArray.Create(1, 1);
                resultHessYX = MultidimensionalArray.Create(1, 1);
                resultHessYY = MultidimensionalArray.Create(1, 1);
            }

            gradientX.Evaluate(jCell, 1, nodeSet, resultGradX);
            gradientY.Evaluate(jCell, 1, nodeSet, resultGradY);

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

        /// <summary>
        /// Algorithm to find an inflection point along a curve in the direction of the gradient
        /// </summary>
        /// <param name="gridData">The corresponding grid</param>
        /// <param name="field">The DG field, which shall be evalauted</param>
        /// <param name="points">
        /// The points along the curve (first point has to be user defined)
        /// Lenghts --> [0]: maxIterations + 1, [1]: 5
        /// [1]: x | y | function values | second derivatives | step sizes
        /// </param>
        /// <param name="converged">Has an inflection point been found?</param>
        /// <param name = "jLocal" >Local cell index of the inflection point (if found)</ param >
        /// < param name="byFlux">Calculation of the first and second order derivatives, default: true</param>
        /// <returns>The amount of iterations needed in order to find the inflection point</returns>
        private int WalkOnCurve(GridData gridData, SinglePhaseField field, MultidimensionalArray points, out bool converged, out int jLocal, bool byFlux = true) {
            // Init
            converged = false;

            // Current (global) point
            double[] currentPoint = new double[] { points[0, 0], points[0, 1] };

            // Compute global cell index of current point
            gridData.LocatePoint(currentPoint, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess);

            // Compute local node set
            NodeSet nodeSet = ShockFindingHelperFunctions.GetLocalNodeSet(gridData, currentPoint, (int)GlobalIndex);

            // Get local cell index of current point
            int j0Grd = gridData.CellPartitioning.i0;
            jLocal = (int)(GlobalIndex - j0Grd);

            // Evaluate the second derivative
            try {
                if (byFlux) {
                    points[0, 3] = SecondDerivativeByFlux(field, jLocal, nodeSet);
                } else {
                    points[0, 3] = ShockFindingHelperFunctions.SecondDerivative(field, jLocal, nodeSet);
                }
            } catch (NotSupportedException) {
                return 0;
            }

            // Evaluate the function
            MultidimensionalArray f = MultidimensionalArray.Create(1, 1);
            field.Evaluate(jLocal, 1, nodeSet, f);
            points[0, 2] = f[0, 0];

            // Set initial step size to 0.5 * h_minGlobal
            // Set the search direction depending on the sign of the curvature
            points[0, 4] = 0.5 * gridData.Cells.h_minGlobal;

            int n = 1;
            while (n < points.Lengths[0]) {
                // Evaluate the gradient of the current point
                double[] gradient;
                if (byFlux) {
                    gradient = NormalizedGradientByFlux(field, jLocal, nodeSet);
                } else {
                    gradient = ShockFindingHelperFunctions.NormalizedGradient(field, jLocal, nodeSet);
                }

                // Compute new point along curve
                currentPoint[0] = currentPoint[0] + gradient[0] * points[n - 1, 4];
                currentPoint[1] = currentPoint[1] + gradient[1] * points[n - 1, 4];

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
                            nodeSet = ShockFindingHelperFunctions.GetLocalNodeSet(gridData, currentPoint, jLocal);
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
                    nodeSet = ShockFindingHelperFunctions.GetLocalNodeSet(gridData, currentPoint, jLocal + j0Grd);
                }

                // Update output
                points[n, 0] = currentPoint[0];
                points[n, 1] = currentPoint[1];

                // Evaluate the function
                f.Clear();
                field.Evaluate(jLocal, 1, nodeSet, f);
                points[n, 2] = f[0, 0];

                // Evaluate the second derivative of the new point
                try {
                    if (byFlux) {
                        points[n, 3] = SecondDerivativeByFlux(field, jLocal, nodeSet);
                    } else {
                        points[n, 3] = ShockFindingHelperFunctions.SecondDerivative(field, jLocal, nodeSet);
                    }
                } catch (NotSupportedException) {
                    break;
                }

                // Check if sign of second derivative has changed
                // Halve the step size and change direction
                if (Math.Sign(points[n, 3]) != Math.Sign(points[n - 1, 3])) {
                    points[n, 4] = -points[n - 1, 4] / 2;
                } else {
                    points[n, 4] = points[n - 1, 4];
                }

                n++;

                // Termination criterion
                // Remark: n has already been incremented!
                double[] diff = new double[currentPoint.Length];
                diff[0] = points[n - 1, 0] - points[n - 2, 0];
                diff[1] = points[n - 1, 1] - points[n - 2, 1];
                if (diff.L2Norm() < 1e-12) {
                    converged = true;
                    break;
                }
            }

            return n;
        }
    }

    /// <summary>
    /// Static methods for locating the inflection points of a DG Field
    /// in two dimensions. Can be used to find the shock location. 
    /// </summary>
    public static class ShockFindingHelperFunctions {
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
        /// Returns all cells with the respective artificial viscosity value
        /// if the value is larger than zero
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="avField"></param>
        /// <returns><see cref="MultidimensionalArray"/>, first index: cellIndex, second index: artificial viscosity value</returns>
        public static MultidimensionalArray GetAVMeanValues(GridData gridData, SinglePhaseField avField) {
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

        //public static MultidimensionalArray SortPoints(GridData gridData, SinglePhaseField field, MultidimensionalArray points, bool byFlux = false) {
        //    double[] sortedXCoords = new double[points.Lengths[0]];
        //    double[] sortedYCoords = new double[points.Lengths[0]];
        //    int keptEntries = 0;

        //    for (int i = 0; i < points.Lengths[0]; i++) {
        //        // Compute global cell index of the point
        //        double[] currentPoint = points.ExtractSubArrayShallow(i, 0, -1).To1DArray();
        //        gridData.LocatePoint(currentPoint, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess);

        //        // Compute local node set
        //        NodeSet nodeSet = GetLocalNodeSet(gridData, currentPoint, (int)GlobalIndex);

        //        // Get local cell index of current point
        //        int j0Grd = gridData.CellPartitioning.i0;
        //        int jLocal = (int)(GlobalIndex - j0Grd);

        //        // Calculate secondDerivative
        //        double secondDerivative;
        //        if (byFlux) {
        //            secondDerivative = SecondDerivativeByFlux(field, jLocal, nodeSet,);
        //        } else {
        //            secondDerivative = SecondDerivative(field, jLocal, nodeSet);
        //        }

        //        // Select points where second derivative is larger than zero
        //        if (secondDerivative > 0.0) {
        //            sortedXCoords[keptEntries] = currentPoint[0];
        //            sortedYCoords[keptEntries] = currentPoint[1];
        //            keptEntries++;
        //        }
        //    }

        //    // Resize final array
        //    MultidimensionalArray result = MultidimensionalArray.Create(keptEntries, points.Lengths[1], points.Lengths[2]);
        //    for (int i = 0; i < keptEntries; i++) {
        //        result[i, 0, 0] = sortedXCoords[i];
        //        result[i, 0, 1] = sortedYCoords[i];
        //    }

        //    return result;
        //}

        /// <summary>
        /// Performs the patch recovery on a given DG field
        /// </summary>
        /// <param name="input">A <see cref="SinglePhaseField"/></param>
        /// <returns>A patch recovered <see cref="SinglePhaseField"/> with a degree of inputDegreee + 2</returns>
        public static SinglePhaseField PatchRecovery(SinglePhaseField input) {
            Console.WriteLine(String.Format("Patch recovery of field {0} started...", input.Identification));

            Basis prcBasis = new Basis(input.GridDat, input.Basis.Degree + 2);
            L2PatchRecovery prc = new L2PatchRecovery(input.Basis, prcBasis, CellMask.GetFullMask(input.GridDat), RestrictToCellMask: true);
            SinglePhaseField prcField = new SinglePhaseField(prcBasis, input.Identification + "_prc");
            prc.Perform(prcField, input);

            Console.WriteLine(String.Format("finished", input.Identification));

            return prcField;
        }

        public static SinglePhaseField ReconstructLevelSetField(SinglePhaseField field, MultidimensionalArray points) {
            // Init
            IGridData gridData = field.GridDat;

            // Evaluate gradient
            Basis basis = new Basis(field.GridDat, field.Basis.Degree);
            SinglePhaseField gradientX = new SinglePhaseField(basis, "gradientX");
            SinglePhaseField gradientY = new SinglePhaseField(basis, "gradientY");
            gradientX.DerivativeByFlux(1.0, field, d: 0);
            gradientY.DerivativeByFlux(1.0, field, d: 1);
            //gradientX = PatchRecovery(gradientX);
            //gradientY = PatchRecovery(gradientY);
            FieldEvaluation ev = new FieldEvaluation((GridData)gridData);
            MultidimensionalArray GradVals = MultidimensionalArray.Create(points.GetLength(0), 2);
            ev.Evaluate(1.0, new DGField[] { gradientX, gradientY }, points, 0.0, GradVals);

            // Level set reconstruction
            Console.WriteLine(String.Format("Reconstruction of level set started...", field.Identification));
            int count = 0;
            Func<double[], double> func = delegate (double[] X) {
                double minDistSigned = double.MaxValue;
                int iMin = int.MaxValue;
                double closestPointOnInterfaceX = double.MaxValue;
                double closestPointOnInterfaceY = double.MaxValue;

                double x = X[0];
                double y = X[1];

                for (int i = 0; i < points.Lengths[0]; i++) {
                    double currentPointX = points[i, 0];
                    double currentPointY = points[i, 1];

                    double deltaX = x - currentPointX;
                    double deltaY = y - currentPointY;

                    double dist = Math.Sqrt(deltaX * deltaX + deltaY * deltaY);

                    if (dist <= minDistSigned) {
                        iMin = i;
                        minDistSigned = dist;
                        closestPointOnInterfaceX = currentPointX;
                        closestPointOnInterfaceY = currentPointY;
                    }
                }

                //// Compute global cell index of quadrature node
                //gridData.LocatePoint(X, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess);

                //// Compute local node set
                //NodeSet nodeSet = GetLocalNodeSet(gridData, X, (int)GlobalIndex);

                //// Evaluate gradient
                //// Get local cell index of quadrature node
                //int j0Grd = gridData.CellPartitioning.i0;
                //int j0 = (int)(GlobalIndex - j0Grd);

                //int Len = 1;
                //MultidimensionalArray resultGradX = MultidimensionalArray.Create(1, 1);
                //MultidimensionalArray resultGradY = MultidimensionalArray.Create(1, 1);
                //gradientX.Evaluate(j0, Len, nodeSet, resultGradX);
                //gradientY.Evaluate(j0, Len, nodeSet, resultGradY);

                int sign = Math.Sign((x - closestPointOnInterfaceX) * GradVals[iMin, 0] + (y - closestPointOnInterfaceY) * GradVals[iMin, 1]);
                minDistSigned *= sign;

                //Console.WriteLine(String.Format("Quadrature point #{0}: ({1}, {2}), interface point: ({3}, {4})", count, X[0], X[1], closestPointOnInterfaceX, closestPointOnInterfaceY));
                count++;

                return minDistSigned;
            };

            // DG level set field
            SinglePhaseField DGLevelSet = new SinglePhaseField(field.Basis, "levelSet_recon");
            DGLevelSet.ProjectField(func.Vectorize());

            Console.WriteLine("finished");
            return DGLevelSet;
        }

        /// <summary>
        /// Creates a continuous version of a DG level set field
        /// </summary>
        /// <param name="DGLevelSet">The input DG level set field</param>
        /// <returns>A continuous level set <see cref="SinglePhaseField"/> with one order higher than the input</returns>
        public static SinglePhaseField ContinuousLevelSet(SinglePhaseField DGLevelSet, MultidimensionalArray points) {
            // Init
            IGridData gridData = DGLevelSet.GridDat;
            SinglePhaseField continuousLevelSet = new SinglePhaseField(new Basis(gridData, DGLevelSet.Basis.Degree + 1), DGLevelSet.Identification + "_cont");

            // Create narrow band
            double[] jCells = points.ExtractSubArrayShallow(-1, 3).To1DArray();
            List<int> allNeighbours = new List<int>();
            foreach (int cell in jCells) {
                gridData.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] cellNeighbours, out int[] connectingEntities);
                allNeighbours.Add(cell);
                allNeighbours.AddRange(cellNeighbours);
            }
            List<int> narrowBandCells = allNeighbours.Distinct().ToList();
            BitArray bitMaskNarrowBand = new BitArray(gridData.iLogicalCells.NoOfLocalUpdatedCells);
            foreach (int cell in narrowBandCells) {
                bitMaskNarrowBand[cell] = true;
            }
            CellMask narrowBand = new CellMask(gridData, bitMaskNarrowBand);

            // Create positive mask
            BitArray bitMaskPos = new BitArray(gridData.iLogicalCells.NoOfLocalUpdatedCells);
            for (int i = 0; i < gridData.iLogicalCells.NoOfLocalUpdatedCells; i++) {
                if (DGLevelSet.GetMeanValue(i) > 0 && !bitMaskNarrowBand[i]) {
                    bitMaskPos[i] = true;
                }
            }
            CellMask posMask = new CellMask(gridData, bitMaskPos);

            // Continuity projection
            ContinuityProjection continuityProjection = new ContinuityProjection(continuousLevelSet.Basis, DGLevelSet.Basis, gridData, ContinuityProjectionOption.ConstrainedDG);
            continuityProjection.MakeContinuous(DGLevelSet, continuousLevelSet, narrowBand, posMask);

            Console.WriteLine(String.Format("Continuity projection of field {0} finished", DGLevelSet.Identification));

            return continuousLevelSet;
        }

        /// <summary>
        /// Reads a text where a clustering is saved
        /// Index 0: x-coordinate, index 1: y-coordinate, index 2: field value, index 3: cluster number
        /// </summary>
        /// <param name="path">Path of the file to read</param>
        public static MultidimensionalArray ReadTextFile(string path) {
            string[] lines = File.ReadAllLines(path);
            int cols = lines[0].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries).Length;
            MultidimensionalArray result = MultidimensionalArray.Create(lines.Length, cols);

            for (int i = 0; i < lines.Length; i++) {
                for (int j = 0; j < cols; j++) {
                    result[i, j] = Convert.ToDouble(lines[i].Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries)[j]);
                }
            }

            return result;
        }

        public static DGField Find(this IEnumerable<DGField> fields, string name) {
            return fields.Where(f => f.Identification == name).SingleOrDefault();
        }
    }

    public enum SeedingSetup {
        av = 0,
        av3x3
    }
}
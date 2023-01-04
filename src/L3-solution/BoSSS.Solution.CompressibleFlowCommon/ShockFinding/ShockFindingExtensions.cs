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
using BoSSS.Foundation.ConstrainedDGprojection;
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
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockFinding {

    public static class ShockFindingExtensions {

        public static MultidimensionalArray LoadResults(string path) {
            MultidimensionalArray[] parts = new MultidimensionalArray[5];
            for (int i = 0; i < 5; i++) {
                parts[i] = IMatrixExtensions.LoadFromTextFile(path + "Results_" + i + ".txt");
            }

            MultidimensionalArray result = MultidimensionalArray.Create(parts[0].Lengths[0], parts[0].Lengths[1], 5);
            for (int i = 0; i < result.Lengths[2]; i++) {
                result.ExtractSubArrayShallow(-1, -1, i).Acc(1.0, parts[i]);
            }
            return result;
        }

        public static MultidimensionalArray LoadResultsExtended(string path) {
            return IMatrixExtensions.LoadFromTextFile(path + "ResultsExtended.txt");
        }

        public static void SaveResults(this MultidimensionalArray input, string path) {
            if (input.Dimension == 2) {
                input.SaveToTextFile(path + "ResultsExtended.txt");
            } else if (input.Dimension == 3) {
                for (int i = 0; i < input.Lengths.Last(); i++) {
                    input.ExtractSubArrayShallow(-1, -1, i).SaveToTextFile(path + "Results_" + i + ".txt");
                }
            } else {
                throw new NotSupportedException("This MdA cannot be saved.");
            }
        }

        public static void EliminateNonConvergedPoints(MultidimensionalArray input, MultidimensionalArray inputExtended, out MultidimensionalArray result, out MultidimensionalArray resultExtended) {
            // input            [0]: x        [1]: y             [2]: f       [3] secondDerivative    [4] stepSize
            // inputExtended    [0]: iter     [1]: converged     [2]: jCell

            Console.WriteLine("EliminateNonConvergedPoints: START");

            int numOfPoints = input.Lengths[0];
            int[] convergedCells = new int[numOfPoints];
            int count = 0;
            for (int i = 0; i < numOfPoints; i++) {
                if (inputExtended[i, 1] > 0.0) {
                    convergedCells[count] = i;
                    count++;
                }
            }
            Array.Resize(ref convergedCells, count);

            result = MultidimensionalArray.Create(count, input.Lengths[1], input.Lengths[2]);
            resultExtended = MultidimensionalArray.Create(count, inputExtended.Lengths[1]);

            for (int i = 0; i < convergedCells.Length; i++) {
                // If converged, copy all info from input to output array
                int cell = convergedCells[i];
                result.ExtractSubArrayShallow(i, -1, -1).Acc(1.0, input.ExtractSubArrayShallow(cell, -1, -1));
                resultExtended.ExtractSubArrayShallow(i, -1).Acc(1.0, inputExtended.ExtractSubArrayShallow(cell, -1));
            }

            Console.WriteLine("EliminateNonConvergedPoints: END");
        }

        public static double[] GetFinalFunctionValues(MultidimensionalArray input, MultidimensionalArray iterationsNeeded) {
            double[] result = new double[input.Lengths[0]];
            for (int i = 0; i < input.Lengths[0]; i++) {
                result[i] = input[i, (int)iterationsNeeded[i] - 1, 2];
            }

            return result;
        }

        public static string[] CreateDirectories(string mainPath, string directoryName, List<ISessionInfo> sessions) {
            string[] sessionPaths = new string[sessions.Count()];
            for (int i = 0; i < sessions.Count(); i++) {
                sessionPaths[i] = mainPath + sessions.ElementAt(i).Name + @"\" + directoryName + @"\";
                if (!Directory.Exists(sessionPaths[i])) {
                    System.IO.Directory.CreateDirectory(sessionPaths[i]);
                }
            }

            return sessionPaths;
        }

        public static void EmptyDirectories(string mainPath, string directoryName, List<ISessionInfo> sessions) {
            for (int i = 0; i < sessions.Count(); i++) {
                string sessionPath = mainPath + sessions.ElementAt(i).Name + @"\" + directoryName + @"\";
                if (!Directory.Exists(sessionPath)) {
                    return;
                } else {
                    DirectoryInfo sessionDirectory = new DirectoryInfo(sessionPath);

                    // Delete all directories in the current directory
                    foreach (DirectoryInfo dir in sessionDirectory.EnumerateDirectories()) {
                        dir.Delete(true);
                    }

                    // Delete all files in the current directory
                    foreach (FileInfo file in sessionDirectory.EnumerateFiles()) {
                        file.Delete();
                    }
                }
            }
        }

        public static NodeSet GetLocalNodeSet(GridData gridData, double[] globalPoint, long globalCellIndex) {
            int D = 2;
            int localCellIndex = gridData.CellPartitioning.Global2Local(globalCellIndex);
            MultidimensionalArray GlobalVerticesIn = MultidimensionalArray.CreateWrapper(globalPoint, 1, D);
            MultidimensionalArray LocalVerticesOut = MultidimensionalArray.CreateWrapper(new double[D], 1, D);

            gridData.TransformGlobal2Local(GlobalVerticesIn, LocalVerticesOut, localCellIndex, NewtonConvergence: null);

            return new NodeSet(gridData.Cells.GetRefElement(localCellIndex), LocalVerticesOut, true);
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
        /// <param name="gridData">The needed <see cref="GridData"/></param>
        /// <param name="avField">The artificial visocsity <see cref="SinglePhaseField"/></param>
        /// <returns> <see cref="MultidimensionalArray"/>
        /// Lenghts --> [0]: number of points (AV > 0), [1]: 2
        /// [1] --> [0]: cellIndex, [2:] AV value
        /// </returns>
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

        /// <summary>
        /// Reconstructs a <see cref="BoSSS.Foundation.XDG.LevelSet"/> from a given set of points on a given DG field
        /// </summary>
        /// <param name="field">The DG field to work with</param>
        /// <param name="points"> <see cref="MultidimensionalArray"/>
        ///  Lengths --> [0]: numOfPoints, [1]: 2
        ///  [1] --> [0]: x, [1]: y
        /// </param>
        public static SinglePhaseField ReconstructLevelSetField(SinglePhaseField field, MultidimensionalArray points) {
            // Init
            IGridData gridData = field.GridDat;

            // Evaluate gradient
            SinglePhaseField gradientX = new SinglePhaseField(field.Basis, "gradientX");
            SinglePhaseField gradientY = new SinglePhaseField(field.Basis, "gradientY");
            gradientX.DerivativeByFlux(1.0, field, d: 0);
            gradientY.DerivativeByFlux(1.0, field, d: 1);
            //gradientX = PatchRecovery(gradientX);
            //gradientY = PatchRecovery(gradientY);
            FieldEvaluation ev = new FieldEvaluation((GridData)gridData);
            MultidimensionalArray GradVals = MultidimensionalArray.Create(points.GetLength(0), 2);
            ev.Evaluate(1.0, new DGField[] { gradientX, gradientY }, points, 0.0, GradVals);

            // Level set reconstruction
            Console.WriteLine("Reconstruction of field levelSet started...");
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
        /// <param name="jCells">Local cell indices</param>      
        /// <returns>A continuous level set <see cref="SinglePhaseField"/> with one order higher than the input</returns>
        public static SinglePhaseField ContinuousLevelSet(SinglePhaseField DGLevelSet, double[] jCells) {
            Console.WriteLine(String.Format("Continuity projection of field {0} started...", DGLevelSet.Identification));

            // Init
            IGridData gridData = DGLevelSet.GridDat;
            SinglePhaseField continuousLevelSet = new SinglePhaseField(new Basis(gridData, DGLevelSet.Basis.Degree + 1), DGLevelSet.Identification + "_cont");

            // Create narrow band
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


            Console.WriteLine("finished");

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
    }
}
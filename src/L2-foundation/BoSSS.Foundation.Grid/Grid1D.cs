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

using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// onedimensional grid
    /// </summary>
    [Serializable]
    public sealed class Grid1D : GridCommons {

        /// <summary>
        /// constructs an empty 1D grid.
        /// </summary>
        public Grid1D()
            : base(new RefElement[] { Line.Instance }, new RefElement[] { Point.Instance }) {
        }
        

        /// <summary>
        /// creates exponentially distributed nodes between two points
        /// </summary>
        /// <param name="loBound">minimum;</param>
        /// <param name="hiBound">maximum;</param>
        /// <param name="N">number of nodes, i.e. length of the returend array</param>
        /// <param name="alpha">
        /// stretching; must be grater than 0; if grater than one distribution of points gets finer
        /// around <paramref name="loBound"/> and if smaller it becomes coarser.
        /// </param>
        /// <returns></returns>
        static public double[] ExponentialSpaceing(double loBound, double hiBound, int N, double alpha) {
            if (loBound >= hiBound)
                throw new ArgumentException("loBound >= hiBound");
            if (alpha < 0.0)
                throw new ArgumentException("stretching out of range.");

            if (Math.Abs(alpha - 1.0) <= 0.0001) {
                return GenericBlas.Linspace(loBound, hiBound, N);
            } else {
                double[] ret = new double[N];
                N--;

                for (int i = 0; i <= N; i++) {
                    ret[i] = loBound + ((Math.Pow(alpha, i) - 1.0) / (Math.Pow(alpha, N) - 1.0)) * (hiBound - loBound);
                }

                return ret;
            }
        }

        /// <summary>
        /// creates points with hyperbolic tangent distribution.
        /// </summary>
        /// <param name="loBound">minimum value of interval</param>
        /// <param name="hiBound">maximum value of interval</param>
        /// <param name="N">number of nodes, i.e. length of the returend array</param>
        /// <param name="alpha">stretching; must be grater than zero</param>
        /// <param name="loBoundStP">true, if finer distribution towards <paramref name="loBound"/> is desired</param>
        ///  
        /// <returns></returns>
        static public double[] TanhSpacing(double loBound, double hiBound, int N, double alpha, Boolean loBoundStP) {
            if (loBound >= hiBound)
                throw new ArgumentException("loBound >= hiBound");
            if (alpha <= 0)
                throw new ArgumentException("stretching out of range.");

            double[] ret = new double[N];
            N--;
            if (loBoundStP)
                for (int i = 0; i <= N; i++)
                    ret[i] = loBound + (1 + Math.Tanh(alpha * ((double)i / N - 1)) / Math.Tanh(alpha)) * (hiBound - loBound);
            else
                for (int i = 0; i <= N; i++)
                    ret[i] = loBound + Math.Tanh(alpha * ((double)i / N)) / Math.Tanh(alpha) * (hiBound - loBound);

            return ret;
        }

        /// <summary>
        /// creates points with double sided symmetric hyperbolic tangent distribution.
        /// </summary>
        /// <param name="loBound">minimum value of interval</param>
        /// <param name="hiBound">maximum value of interval</param>
        /// <param name="N">number of nodes, i.e. length of the returend array</param>
        /// <param name="alpha">stretching; must be grater than zero</param>
        /// <returns></returns>
        /// <remarks>finer distribution near the two ends</remarks>
        static public double[] TanhSpacing_DoubleSided(double loBound, double hiBound, int N, double alpha) {
            if (loBound >= hiBound)
                throw new ArgumentException("loBound >= hiBound");
            if (alpha <= 0)
                throw new ArgumentException("stretching out of range.");

            double[] ret = new double[N];
            N--;
            for (int i = 0; i <= N; i++)
                ret[i] = loBound + 0.5 * (1 + Math.Tanh(alpha * ((double)i / N - 0.5)) / Math.Tanh(alpha / 2)) * (hiBound - loBound);
            return ret;
        }

        /// <summary>
        /// constructs a 1D-grid.
        /// </summary>
        /// <param name="Nodes">
        /// nodes of the grid.
        /// </param>
        /// <param name="periodic"></param>
        static public Grid1D LineGrid(double[] Nodes, bool periodic = false) {
            return LineGrid(new double[][] { Nodes }, periodic);
        }

        /// <summary>
        /// constructs a 1D-grid consisting of multiple (not connected) segments
        /// </summary>
        /// <param name="nodeSets">
        /// a collection of node segments;
        /// if running with more than one MPI process, only the input on rank 0 is considered
        /// </param>
        /// <param name="periodic"></param>
        static public Grid1D LineGrid(ICollection<double[]> nodeSets, bool periodic = false) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                
                int size, rank;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);

                // Don't forget minus one (sets contain _nodes_, not cells)!
                int globalNumberOfCells = nodeSets.Select(s => s.Length - 1).Sum();

                int globalIndexFirstLocalCell = globalNumberOfCells * rank / size;
                int globalIndexFirstExternalCell = globalNumberOfCells * (rank + 1) / size;

                Grid1D grid = new Grid1D();

                if (Math.Abs(globalIndexFirstLocalCell - globalIndexFirstExternalCell) <= 0) {
                    // No cells on current rank (too few cells?)
                    return grid;
                }

                int localNoOfCells = globalIndexFirstExternalCell - globalIndexFirstLocalCell;
                grid.Cells = new Cell[localNoOfCells];

                byte periodicTag = 0;
                if (periodic) {
                    double[][] left = { new double[] { nodeSets.First().First() } };
                    double[][] right = { new double[] { nodeSets.Last().Last() } };
                    grid.ConstructPeriodicEdgeTrafo(
                        right, new double[] { 1.0 }, left, new double[] { 1.0 }, out periodicTag);
                    grid.EdgeTagNames.Add(periodicTag, "Periodic-X");
                }
                
                int globalCellIndex = -1;
                int localCellIndex = 0;
                foreach (double[] nodesInSegment in nodeSets) {
                    int noOfCellsInSegment = nodesInSegment.Length - 1;

                    for (int i = 0; i < noOfCellsInSegment; i++) {
                        globalCellIndex++;

                        // Handled on other processor
                        if (globalCellIndex < globalIndexFirstLocalCell) {
                            continue;
                        } else if (globalCellIndex >= globalIndexFirstExternalCell) {
                            break;
                        }

                        if (nodesInSegment[i] >= nodesInSegment[i + 1]) {
                            throw new ArgumentException(
                                "Nodes must be provided in strictly ascending order");
                        }

                        Cell cell = new Cell();
                        cell.GlobalID = globalCellIndex;
                        cell.Type = CellType.Line_2;

                        cell.TransformationParams = MultidimensionalArray.Create(2, 1);
                        cell.TransformationParams[0, 0] = nodesInSegment[i];
                        cell.TransformationParams[1, 0] = nodesInSegment[i + 1];

                        cell.NodeIndices = new int[] {
                            globalCellIndex,
                            globalCellIndex + 1
                        };

                        if (periodic) {
                            if (globalCellIndex == 0) {
                                (new CellFaceTag() {
                                    EdgeTag = periodicTag,
                                    PeriodicInverse = true,
                                    ConformalNeighborship = true,
                                    FaceIndex = (int)Line.Edge.Left,
                                    NeighCell_GlobalID = globalNumberOfCells - 1
                                }).AddToArray(ref cell.CellFaceTags);
                            } else if (globalCellIndex == globalNumberOfCells - 1) {
                                (new CellFaceTag() {
                                    EdgeTag = periodicTag,
                                    PeriodicInverse = false,
                                    ConformalNeighborship = true,
                                    FaceIndex = (int)Line.Edge.Right,
                                    NeighCell_GlobalID = 0
                                }).AddToArray(ref cell.CellFaceTags);
                            }
                        }

                        grid.Cells[localCellIndex] = cell;

                        localCellIndex++;
                    }
                }

                return grid;
            }
        }
    }
}

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
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockFinding {

    public class LevelSetReconstruction {

        private readonly string sessionPath;
        private readonly ISessionInfo session;
        private readonly MultidimensionalArray input;
        private readonly MultidimensionalArray inputExtended;

        private int clusteringCount = 0;

        private readonly List<MultidimensionalArray> _clusterings = new List<MultidimensionalArray>();
        /// <summary>
        /// List of all clusterings in the order of their creation
        /// </summary>
        public List<MultidimensionalArray> Clusterings {
            get {
                if (_clusterings == null) {
                    throw new NotImplementedException("List of clusterings is empty.");
                }
                return _clusterings;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="sessionPath">Path where everything is stored</param>
        /// <param name="session">The session from the database</param>
        /// <param name="input">
        /// Lenghts --> [0]: numOfPoints, [1]: maxIterations + 1, [2]: 5
        /// [2]: x | y | function values | second derivatives | step sizes
        /// </param>
        /// <param name="inputExtended">
        /// Lenghts --> [0]: numOfPoints, [1]: 3
        /// [1]: IterationsNeeded | Converged | jCell
        /// </param>
        public LevelSetReconstruction(string sessionPath, ISessionInfo session, MultidimensionalArray input, MultidimensionalArray inputExtended) {
            this.sessionPath = sessionPath;
            this.session = session;
            this.input = input;
            this.inputExtended = inputExtended;
        }

        /// <summary>
        /// Create clustering based on the density
        /// </summary>
        /// <param name="numOfClusters">Needed by <see cref="Kmeans"/></param>
        /// <param name="initialMeans">Needed by <see cref="Kmeans"/></param>
        /// <returns>
        /// Clustering as <see cref="MultidimensionalArray"/>
        /// [0]: x, [1]: y, [2]: data, [3]: cellToCluster (e.g. cell 0 is in cluster 1), [4]: local cell index
        /// </returns>
        public MultidimensionalArray CreateClustering_Density(int numOfClusters, double[] initialMeans) {
            Console.WriteLine("CreateClustering_Density: START");

            double[] data = ShockFindingExtensions.GetFinalFunctionValues(input, inputExtended.ExtractSubArrayShallow(-1, 0));
            Kmeans kmeans = new Kmeans(data, numOfClusters, initialMeans);
            int[] cellToCluster = kmeans.Cluster();

            MultidimensionalArray clustering = MultidimensionalArray.Create(data.Length, 5);
            for (int i = 0; i < input.Lengths[0]; i++) {
                clustering[i, 0] = input[i, (int)inputExtended[i, 0] - 1, 0];      // x
                clustering[i, 1] = input[i, (int)inputExtended[i, 0] - 1, 1];      // y
                clustering[i, 2] = data[i];                                        // data value
                clustering[i, 3] = cellToCluster[i];                               // cellToCluster (e.g. cell 0 is in cluster 1)
                clustering[i, 4] = inputExtended[i, 2];                            // local cell index
            }
            _clusterings.Add(clustering);

            Console.WriteLine("CreateClustering_Density: END");

            return clustering;
        }

        /// <summary>
        /// Create clustering based on the artificial viscosity (mean values)
        /// </summary>
        /// <param name="inputClustering">Input data which has to be a previous clustering</param>
        /// <param name="numOfClusters">Needed by <see cref="Kmeans"/></param>
        /// <param name="initialMeans">Needed by <see cref="Kmeans"/></param>
        /// <returns>
        /// Clustering as <see cref="MultidimensionalArray"/>
        /// [0]: x, [1]: y, [2]: data, [3]: cellToCluster (e.g. cell 0 is in cluster 1), [4]: local cell index
        /// </returns>
        public MultidimensionalArray CreateClustering_AV(MultidimensionalArray inputClustering, int numOfClusters, double[] initialMeans) {
            Console.WriteLine("CreateClustering_AV: START");

            // Get AV values
            var avField = this.session.Timesteps.Last().Fields.Where(f => f.Identification == "artificialViscosity").SingleOrDefault();
            int numOfPoints = inputClustering.Lengths[0];
            double[] data = new double[numOfPoints];
            for (int i = 0; i < data.Length; i++) {
                data[i] = avField.GetMeanValue((int)inputClustering[i, 4]);
            }

            // Kmeans
            Kmeans kmeans = new Kmeans(data, numOfClusters, initialMeans);
            int[] cellToCluster = kmeans.Cluster();

            // Store values
            MultidimensionalArray clustering = MultidimensionalArray.Create(data.Length, inputClustering.Lengths[1]);
            for (int i = 0; i < numOfPoints; i++) {
                clustering[i, 0] = inputClustering[i, 0];      // x
                clustering[i, 1] = inputClustering[i, 1];      // y
                clustering[i, 2] = data[i];                    // data value
                clustering[i, 3] = cellToCluster[i];           // cellToCluster (e.g. cell 0 is in cluster 1)
                clustering[i, 4] = inputClustering[i, 4];      // local cell index 
            }
            _clusterings.Add(clustering);

            Console.WriteLine("CreateClustering_AV: END");

            return clustering;
        }

        /// <summary>
        /// Create clustering without boundary cells
        /// </summary>
        /// <param name="inputClustering">Input data which has to be a previous clustering</param>
        /// <returns>
        /// Clustering as <see cref="MultidimensionalArray"/>
        /// [0]: x, [1]: y, [2]: data, [3]: cellToCluster (e.g. cell 0 is in cluster 1), [4]: local cell index
        public MultidimensionalArray CreateClustering_Boundary(MultidimensionalArray inputClustering) {
            Console.WriteLine("CreateClustering_Boundary: START");

            var gridData = (GridData)session.Timesteps.Last().Fields.First().GridDat;
            BitArray isBoundaryCell = gridData.GetBoundaryCells().GetBitMask();

            // Store values
            int numOfPoints = inputClustering.Lengths[0];
            int[] internalCells = new int[numOfPoints];
            int count = 0;
            for (int i = 0; i < numOfPoints; i++) {
                if (!isBoundaryCell[(int)inputClustering[i, 4]]) {
                    internalCells[count] = i;
                    count++;
                }
            }
            Array.Resize(ref internalCells, count);

            MultidimensionalArray clustering = MultidimensionalArray.Create(count, inputClustering.Lengths[1]);
            for (int i = 0; i < internalCells.Length; i++) {
                int cell = internalCells[i];
                clustering.ExtractSubArrayShallow(i, -1).Acc(1.0, inputClustering.ExtractSubArrayShallow(cell, -1));
            }
            _clusterings.Add(clustering);

            Console.WriteLine("CreateClustering_Boundary: END");

            return clustering;
        }

        /// <summary>
        /// Select a specific cluster from a clustering and extract all data
        /// </summary>
        /// <param name="clustering">The underlying clustering</param>
        /// <param name="clusterToSelect">The cluster to extract</param>
        /// <returns>New clustering which contains only the data from the selected cluster</returns>
        public MultidimensionalArray SelectCluster(MultidimensionalArray clustering, int clusterToSelect) {
            // clustering.Lengths -->  [0]: numOfPoints      [1]: 4
            // clustering[1] -->       [0]: x                [1]: y          [2]: function value       [3]: cellToCluster       [4]: local cell index

            int numOfPoints = clustering.Lengths[0];
            int[] cellsInCluster = new int[numOfPoints];
            int count = 0;
            for (int i = 0; i < numOfPoints; i++) {
                if (Math.Abs(clustering[i, 3] - clusterToSelect) < 0.1) {
                    cellsInCluster[count] = i;
                    count++;
                }
            }
            Array.Resize(ref cellsInCluster, count);

            MultidimensionalArray result = MultidimensionalArray.Create(count, clustering.Lengths[1]);
            for (int i = 0; i < cellsInCluster.Length; i++) {
                int cell = cellsInCluster[i];
                result.ExtractSubArrayShallow(i, -1).Acc(1.0, clustering.ExtractSubArrayShallow(cell, -1));
            }

            return result;
        }

        /// <summary>
        /// Save a clustering to a text file
        /// </summary>
        /// <param name="clustering">[0]: x, [1]: y, [2]: data, [3]: cellToCluster (e.g. cell 0 is in cluster 1), [4]: local cell index</param>
        /// <param name="path">Optional path</param>
        public void SaveClusteringToTextFile(MultidimensionalArray clustering, string path = null) {
            if (path == null) {
                path = sessionPath;
            }

            using (System.IO.StreamWriter sw = new System.IO.StreamWriter(path + "clustering_" + clusteringCount + ".txt")) {
                string resultLine;
                for (int i = 0; i < clustering.Lengths[0]; i++) {
                    resultLine = clustering[i, 0]
                        + "\t" + clustering[i, 1]
                        + "\t" + clustering[i, 2]
                        + "\t" + clustering[i, 3];
                    sw.WriteLine(resultLine);
                }
                sw.Flush();
            }

            clusteringCount++;
        }
    }
}
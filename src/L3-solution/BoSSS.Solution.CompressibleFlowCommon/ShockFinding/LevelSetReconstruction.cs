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

        private readonly List<MultidimensionalArray> _clusterings = new List<MultidimensionalArray>();
        public List<MultidimensionalArray> Clusterings {
            get {
                if (_clusterings == null) {
                    throw new NotImplementedException("List of clusterings is empty.");
                }
                return _clusterings;
            }
        }

        public LevelSetReconstruction(string sessionPath, ISessionInfo session, MultidimensionalArray input, MultidimensionalArray inputExtended) {
            this.sessionPath = sessionPath;
            this.session = session;
            this.input = input;
            this.inputExtended = inputExtended;
        }

        private double[] SortOutNonConverged(double[] data) {
            double[] result = new double[data.Length];

            int count = 0;
            for (int i = 0; i < input.Lengths[0]; i++) {
                if (inputExtended[i, 1] > 0) {
                    result[count] = data[i];
                    count++;
                }
            }

            Array.Resize(ref result, count);
            return result;
        }

        public MultidimensionalArray CreateClustering(int numOfClusters, double[] initialMeans, double[] data = null, bool sortOutNonConverged = false) {
            Console.WriteLine("CLUSTERING: START");

            // Test
            //MultidimensionalArray input = MultidimensionalArray.Create(2, 1, 5);
            //input[0, 0, 0] = 0;
            //input[0, 0, 1] = 0.1;
            //input[0, 0, 2] = 0.2;
            //input[0, 0, 3] = 0.3;
            //input[0, 0, 4] = 0.4;
            //input[1, 0, 0] = 1.0;
            //input[1, 0, 1] = 1.1;
            //input[1, 0, 2] = 1.2;
            //input[1, 0, 3] = 1.3;
            //input[1, 0, 4] = 1.4;

            //MultidimensionalArray inputExtended = MultidimensionalArray.Create(2, 3);
            //inputExtended[0, 0] = 1;
            //inputExtended[0, 1] = 1;
            //inputExtended[0, 2] = 0;
            //inputExtended[1, 0] = 1;
            //inputExtended[1, 1] = -1;
            //inputExtended[1, 2] = 1;

            if (data == null) {
                data = ShockFindingExtensions.GetFinalFunctionValues(input, inputExtended.ExtractSubArrayShallow(-1, 0));
            }
            Kmeans kmeans = new Kmeans(data, numOfClusters, initialMeans);
            int[] cellToCluster = kmeans.Cluster();

            MultidimensionalArray clustering = MultidimensionalArray.Create(data.Length, 4);
            for (int i = 0; i < input.Lengths[0]; i++) {
                clustering[i, 0] = input[i, (int)inputExtended[i, 0] - 1, 0];      // x
                clustering[i, 1] = input[i, (int)inputExtended[i, 0] - 1, 1];      // y
                clustering[i, 2] = data[i];                                        // function value
                clustering[i, 3] = cellToCluster[i];                               // cellToCluster (e.g. cell 0 is in cluster 1)
            }
            _clusterings.Add(clustering);

            Console.WriteLine("CLUSTERING: END");
            return clustering;
        }

        public void SaveClusteringToTextFile(int whichClustering, string path = null) {
            if (path == null) {
                path = sessionPath;
            }

            MultidimensionalArray clustering = _clusterings[whichClustering];
            using (System.IO.StreamWriter sw = new System.IO.StreamWriter(path + "clustering_" + whichClustering + ".txt")) {
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
        }
    }
}
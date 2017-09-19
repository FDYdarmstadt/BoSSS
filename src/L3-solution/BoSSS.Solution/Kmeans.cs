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

using System;
using ilPSP.Utils;
using MPI.Wrappers;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Implementation of the k-means clustering algorithm
    /// </summary>
    class Kmeans {
        private double[] data;
        private double[] oldLocalMeans;
        private int numOfCluster;
        private int maxIerations;
        private int numOfProcs;

        // Getter & Setter
        public int[] Clustered {
            get;
            private set;
        }

        public int[] ClusterCount {
            get;
            private set;
        }

        public double[] Means {
            get;
            private set;
        }

        /// <summary>
        /// Constructor of the k-means algorithm
        /// </summary>
        /// <param name="data">Array to cluster</param>
        /// <param name="NumOfCluster">Number of cluster</param>
        /// <param name="means">Starting mean values for each cluster</param>
        /// <param name="MaxIterations">optional Parameter for the maximum iterations</param>
        /// <remarks>Works only for one dimensional data entries, e.g., double[] not double[][]</remarks>
        public Kmeans(double[] data, int NumOfCluster, double[] means, int MaxIterations=100) {
            this.data = data;
            this.numOfCluster = NumOfCluster;
            this.maxIerations = MaxIterations;
            this.Means = means;
            oldLocalMeans = means;
            if (means.Length != NumOfCluster) throw new ArgumentException("K-means clustering not possible: Number of cluster and number of mean values are not equal");

            // Initialize clustering with -1
            Clustered = new int[data.Length];
            ArrayTools.SetAll(Clustered, -1);

            ClusterCount = new int[NumOfCluster];

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out numOfProcs);
        }

        /// <summary>
        /// Main function: Performs the clustering
        /// </summary>
        /// <returns>Array which cluster each data entry belongs to</returns>
        public int[] Cluster() {
            bool changed = true;
            bool changedMean = true;
            int iter = 0;

            do{
                changed = AssignData(); //Involves MPI
                changedMean = UpdateMeans(); //Involves MPI, to ensure that all processors have the same mean-values
                iter++;
            } while (changed && changedMean && iter <= maxIerations);

            
            Console.WriteLine("K-means finished after "+iter+" iterations");
            return Clustered;
        }

        /// <summary>
        /// Assigns the data to the cluster
        /// </summary>
        /// <returns>True, if at least one data-entry changed its cluster</returns>
        private bool AssignData() {
            bool change = false;
            double[] distance = new double[numOfCluster];
            for (int i = 0; i < data.Length; i++) {
                if (!double.IsNaN(data[i])) { //Occurs for void cells in an IBM simulation
                    for (int j = 0; j < numOfCluster; j++) {
                        distance[j] = Math.Abs(data[i] - Means[j]);
                    }
                    int ClusterIndex = GetClusterIndex(distance);

                    if (ClusterIndex != Clustered[i]) {
                        Clustered[i] = ClusterIndex;
                        change = true;
                    }
                }

            }

            // Global update of change
            if (numOfProcs > 1) {
                int sndBool = 0;
                int rcvBool = 0;
                sndBool = change ? 1 : 0;
                unsafe {
                    csMPI.Raw.Allreduce((IntPtr)(&sndBool), (IntPtr)(&rcvBool), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                }

                change = rcvBool == 0 ? false : true; 
            }

            return change;
        }

        /// <summary>
        /// Returns the id of the cluster with the shortest distance 
        /// </summary>
        /// <param name="distance">Given distance for an data entry to all cluster centers</param>
        /// <returns>Id of the cluster with shortest distance</returns>
        private int GetClusterIndex(double[] distance) {
            int Index=0;
            double MinDistance = distance[0];
            for(int i=1; i< distance.Length; i++){
                if (distance[i] < MinDistance) {
                    MinDistance = distance[i];
                    Index=i;
                }
            }
            return Index;
        }

        /// <summary>
        /// Compute the new mean value of every cluster
        /// </summary>
        /// <returns>True, if at least one mean value has changed</returns>
        private bool UpdateMeans() {
            bool changed = false;
            // Local update
            for (int i = 0; i < numOfCluster; i++) {
                double newMean = 0.0;
                ClusterCount[i]= 0;
                for (int j = 0; j < data.Length; j++) {
                    if (Clustered[j] == i) {
                        newMean += data[j];
                        ClusterCount[i]++;
                    }
                }
                if (ClusterCount[i] != 0) {
                    double mean = newMean / ClusterCount[i];
                    if (Math.Abs(Means[i] - mean) > 1e-10){
                        Means[i] = mean;
                        changed = true;
                    }
                }
            }
            // Global update
            if (numOfProcs > 1) { // MPI only needed for more than one processor

                double[] globalMeans = new double[numOfCluster];
                unsafe {
                    // send = means[]
                    // receive = globalMeans[]
                    fixed (double* pSend = &Means[0], pRcv = &globalMeans[0]) {
                        csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), numOfCluster, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    }
                }

                for (int i = 0; i < numOfCluster; i++) {
                    double mean = globalMeans[i] / numOfProcs;
                    if (Math.Abs(Means[i] - mean) > 1e-10) {
                        changed = true;
                        Means[i] = mean;
                    }
                }
                int sndBool = 0;
                int rcvBool = 0;
                sndBool = changed ? 1 : 0;
                unsafe {
                    csMPI.Raw.Allreduce((IntPtr)(&sndBool), (IntPtr)(&rcvBool), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                }

                changed = rcvBool == 0 ? false : true; 
            }
            return changed;
        }
    }
}
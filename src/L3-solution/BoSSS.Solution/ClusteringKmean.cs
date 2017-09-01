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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.Utils;
using MPI.Wrappers;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Implementation of the k-mean clustering algorithm
    /// </summary>
    class ClusteringKmean {
        private double[] data;
        private double[] oldLocalMeans;
        private int NumOfCluster;
        private int MaxIerations;
        private int NumOfProc;

        // Getter & Setter
        public int[] clusterd {
            get;
            private set;
        }

        public int[] ClusterCount {
            get;
            private set;
        }

        public double[] means {
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
        /// <remarks>Works only for one dimensional data entries, e.g double[] not double[][]</remarks>
        public ClusteringKmean(double[] data, int NumOfCluster, double[] means, int MaxIterations=100) {
            this.data = data;
            this.NumOfCluster = NumOfCluster;
            this.MaxIerations = MaxIterations;
            this.means = means;
            oldLocalMeans = means;
            if (means.Length != NumOfCluster) throw new ArgumentException("K-Mean clustering not possible: Number of cluster and number of Mean-Values are not equal");

            // Initialize clustering with -1
            clusterd = new int[data.Length];
            ArrayTools.SetAll(clusterd, -1);

            ClusterCount = new int[NumOfCluster];

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out NumOfProc);
        }

        /// <summary>
        /// Main function: Performs the clustering
        /// </summary>
        /// <returns>Array to which cluster each data entry belongs</returns>
        public int[] Cluster() {
            bool changed = true;
            bool changedMean = true;
            int iter = 0;

            do{
                changed = AssignData(); //Involves MPI
                changedMean = UpdateMeans(); //Involves MPI, to ensure that all processors have the same mean-values
                iter++;
            } while (changed && changedMean && iter <= MaxIerations);

            
            Console.WriteLine("K-means finished after "+iter+" iterations");
            return clusterd;
        }


        // Helper functions

        /// <summary>
        /// Assigns the data to the cluster
        /// </summary>
        /// <returns>true, if at least one data-entry changed its cluster</returns>
        private bool AssignData() {
            bool change = false;
            double[] distance = new double[NumOfCluster];
            for (int i = 0; i < data.Length; i++) {
                if (!double.IsNaN(data[i])) { //Occurs for void cells in an IBM simulation
                    for (int j = 0; j < NumOfCluster; j++) {
                        distance[j] = Math.Abs(data[i] - means[j]);
                    }
                    int ClusterIndex = GetClusterIndex(distance);

                    if (ClusterIndex != clusterd[i]) {
                        clusterd[i] = ClusterIndex;
                        change = true;
                    }
                }

            }


            // global update of change
            if (NumOfProc > 1) {
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
        /// <returns>true, if at least one mean value changed</returns>
        private bool UpdateMeans() {
            bool changed = false;
            // local update
            for (int i = 0; i < NumOfCluster; i++) {
                double newMean = 0.0;
                ClusterCount[i]= 0;
                for (int j = 0; j < data.Length; j++) {
                    if (clusterd[j] == i) {
                        newMean += data[j];
                        ClusterCount[i]++;
                    }
                }
                if (ClusterCount[i] != 0) {
                    double mean = newMean / ClusterCount[i];
                    if (Math.Abs(means[i] - mean) > 1e-10){
                        means[i] = mean;
                        changed = true;
                    }
                }
            }
            // global update
            if (NumOfProc > 1) { // MPI only needed for more than one processor

                double[] globalMeans = new double[NumOfCluster];
                unsafe {
                    // send = means[]
                    // receive = globalMeans[]
                    fixed (double* pSend = &means[0], pRcv = &globalMeans[0]) {
                        csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), NumOfCluster, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    }
                }

                for (int i = 0; i < NumOfCluster; i++) {
                    double mean = globalMeans[i] / NumOfProc;
                    if (Math.Abs(means[i] - mean) > 1e-10) {
                        changed = true;
                        means[i] = mean;
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

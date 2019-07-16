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

using BoSSS.Foundation.Grid;
using BoSSS.Solution.Timestepping;
using System;

namespace CNS.LoadBalancing {

    /// <summary>
    /// All cells belonging to the same LTS cluster (clusters 0..n) are
    /// assigned to the same performance class. The cluster performing the
    /// largest time-steps corresponds to performance class 0, the cluster with
    /// the smallest time-steps to performance class n.
    /// </summary>
    public class LTSCellClassifier : ICellClassifier {

        /// <summary>
        /// All cells belonging to the same LTS cluster (clusters 0..n) are
        /// assigned to the same performance class. The cluster performing the
        /// largest time-steps corresponds to performance class 0, the cluster with
        /// the smallest time-steps to performance class n.
        /// </summary>
        /// <param name="program"></param>
        /// <returns></returns>
        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            AdamsBashforthLTS ltsTimeStepper = program.TimeStepper as AdamsBashforthLTS;
            if (ltsTimeStepper == null) {
                throw new Exception("LTS cell classifier is only sensible for LTS runs.");
            }

            int noOfClasses = ltsTimeStepper.CurrentClustering.NumberOfClusters;

            int J = program.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellToPerformanceClassMap = new int[J];
            for (int i = 0; i < ltsTimeStepper.CurrentClustering.NumberOfClusters; i++) {
                foreach (Chunk chunk in ltsTimeStepper.CurrentClustering.Clusters[i].VolumeMask) {
                    foreach (int cell in chunk.Elements) {
                        cellToPerformanceClassMap[cell] = i;
                    }
                }
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}

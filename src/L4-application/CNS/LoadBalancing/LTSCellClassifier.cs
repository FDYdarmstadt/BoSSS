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

            int[] cellToPerformanceClassMap = new int[program.Grid.NoOfUpdateCells];
            for (int i = 0; i < ltsTimeStepper.CurrentClustering.NumberOfClusters; i++) {
                int noOfTimesteps = ltsTimeStepper.NumberOfLocalTimeSteps[i];
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

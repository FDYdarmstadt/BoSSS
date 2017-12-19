using CNS;
using CNS.LoadBalancing;
using System;

namespace CNS_MPITests {

    /// <summary>
    /// Random classification of cells for testing purposes
    /// </summary>
    public class RandomCellClassifier : ICellClassifier {

        private int noOfClasses;

        public RandomCellClassifier(int noOfClasses) {
            this.noOfClasses = noOfClasses;
        }

        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            // Use seed to make runs deterministic, without having the same
            // random numbers for all cores
            Random rand = new Random(program.ResLogger.TimeStep + program.MPIRank);

            int[] cellToPerformanceClassMap = new int[program.GridData.Cells.NoOfLocalUpdatedCells];
            for (int i = 0; i < cellToPerformanceClassMap.Length; i++) {
                cellToPerformanceClassMap[i] = rand.Next(0, noOfClasses);
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}

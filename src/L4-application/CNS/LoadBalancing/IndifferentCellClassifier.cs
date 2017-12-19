using System;

namespace CNS.LoadBalancing {

    /// <summary>
    /// All cells have the same performance class (i.e., 0)
    /// </summary>
    [Serializable]
    public class IndifferentCellClassifier : ICellClassifier {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="program"></param>
        /// <returns></returns>
        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            return (1, new int[program.Grid.NoOfUpdateCells]);
        }
    }
}

namespace CNS.LoadBalancing {

    /// <summary>
    /// Classifies cells by assigning them a performance class (a number
    /// greater than or equa to zero)
    /// </summary>
    public interface ICellClassifier {

        /// <summary>
        /// Sorts each cell updated by the local process into a performance
        /// class
        /// </summary>
        /// <param name="program"></param>
        /// <returns>
        /// The total number of performance classes (across all processes) and,
        /// for each cell, its assigned performance class
        /// </returns>
        (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program);
    }
}

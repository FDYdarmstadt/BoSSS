using BoSSS.Foundation.Grid;
using BoSSS.Solution.Timestepping;
using CNS.EquationSystem;
using CNS.ShockCapturing;
using System;

namespace CNS.LoadBalancing {

    /// <summary>
    /// AV cells are 1, all others are 0
    /// </summary>
    public class ArtificialViscosityCellClassifier : ICellClassifier {

        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            if (!program.Control.ActiveOperators.HasFlag(Operators.ArtificialViscosity)) {
                throw new Exception(
                    "The selected cell classifier is only sensible for runs with artificial viscosity");
            }
            
            int noOfClasses = 2;
            int[] cellToPerformanceClassMap = new int[program.Grid.NoOfUpdateCells];

            foreach (Chunk chunk in program.Control.ArtificialViscosityLaw.GetShockedCellMask(program.GridData)) {
                foreach (int cell in chunk.Elements) {
                    cellToPerformanceClassMap[cell] = 1;
                }
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}

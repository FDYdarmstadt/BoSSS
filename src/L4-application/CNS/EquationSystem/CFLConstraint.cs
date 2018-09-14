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
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution;
using ilPSP;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;

namespace CNS.EquationSystem {

    /// <summary>
    /// Abstract base-class for standard implementations of CFL time-step
    /// constraints for explicit time integrators
    /// </summary>
    public abstract class CFLConstraint : TimeStepConstraint {

        /// <summary>
        /// CNS variable fields
        /// </summary>
        protected readonly CNSFieldSet workingSet;


        /// <summary>
        /// Cache for nodal density values
        /// <list type="bullet">
        ///     <item>1st index: Cell index - offset</item>
        ///     <item>2nd index: Node index</item>
        /// </list>
        /// </summary>
        protected MultidimensionalArray densityValues;

        /// <summary>
        /// Caches for nodal momentum values in each coordinate direction. For
        /// each direction:
        /// <list type="bullet">
        ///     <item>1st index: Cell index - offset</item>
        ///     <item>2nd index: Node index</item>
        /// </list>
        /// </summary>
        protected MultidimensionalArray[] momentumValues;

        /// <summary>
        /// Cache for nodal energy values
        /// <list type="bullet">
        ///     <item>1st index: Cell index - offset</item>
        ///     <item>2nd index: Node index</item>
        /// </list>
        /// </summary>
        protected MultidimensionalArray energyValues;

        /// <summary>
        /// Constructs a new constraint
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        public CFLConstraint(IGridData gridData, CNSFieldSet workingSet)
            : base(gridData, workingSet.config.dtMin, workingSet.config.dtMax, workingSet.config.CFLFraction, workingSet.config.Endtime) {
            this.workingSet = workingSet;
        }


        /// <summary>
        /// Determines the maximum admissible step-size in the local chunk of
        /// cells from <paramref name="i0"/> to
        /// <paramref name="i0"/> + <paramref name="Length"/>.
        /// </summary>
        /// <param name="i0">First cell index</param>
        /// <param name="Length">Number of cells</param>
        /// <returns>
        /// The maximum admissible step-size for explicit time integrators in
        /// the local chunk of cells from <paramref name="i0"/> to
        /// <paramref name="i0"/> + <paramref name="Length"/>.
        /// </returns>
        public override double GetLocalStepSize(int i0, int Length) {
            NodeSet evaluationPoints = EvaluationPoints[this.gridData.iGeomCells.GetRefElementIndex(i0)];

            // If called via GetGloballyAdmissibleStepSize, these buffers have
            // already been created with an optimal size
            if (densityValues == null || densityValues.GetLength(0) < Length) {
                int maxNoOfNodesPerCell = EvaluationPoints.Max(c => c.NoOfNodes);

                densityValues = MultidimensionalArray.Create(Length, maxNoOfNodesPerCell);
                momentumValues = new MultidimensionalArray[CNSEnvironment.NumberOfDimensions];
                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    momentumValues[d] = MultidimensionalArray.Create(Length, maxNoOfNodesPerCell);
                }
                energyValues = MultidimensionalArray.Create(Length, maxNoOfNodesPerCell);
            }

            workingSet.Density.Evaluate(i0, Length, evaluationPoints, densityValues, 0.0);
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                workingSet.Momentum[d].Evaluate(i0, Length, evaluationPoints, momentumValues[d], 0, 0.0);
            }
            workingSet.Energy.Evaluate(i0, Length, evaluationPoints, energyValues, 0.0);

            return GetCFLStepSize(i0, Length);
        }

        /// <summary>
        /// Implement this method by determining the maximum admissible
        /// step-size in the local chunk of cells from <paramref name="i0"/> to
        /// <paramref name="i0"/> + <paramref name="Length"/>.
        /// </summary>
        /// <param name="i0">First cell index</param>
        /// <param name="Length">Number of cells</param>
        /// <returns>
        /// The maximum admissible step-size for explicit time integrators in
        /// the local chunk of cells from <paramref name="i0"/> to
        /// <paramref name="i0"/> + <paramref name="Length"/>.
        /// </returns>
        protected abstract double GetCFLStepSize(int i0, int Length);

    }
}

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

using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CutCellQuadrature.TestCases {

    /// <summary>
    /// Generic test case providing some standard behavior.
    /// </summary>
    /// <typeparam name="ShiftType">
    /// The type of the shifts that are acceptable.
    /// </typeparam>
    public abstract class TestCase<ShiftType> : ITestCase where ShiftType : Shift {

        /// <summary>
        /// The index of the current shift <see cref="AllShifts"/>
        /// </summary>
        private int shiftIndex = -1;

        /// <summary>
        /// <see cref="ITestCase.GridSize"/>
        /// </summary>
        public GridSizes GridSize {
            get;
            private set;
        }

        /// <summary>
        /// <see cref="ITestCase.GridType"/>
        /// </summary>
        public GridTypes GridType {
            get;
            private set;
        }

        /// <summary>
        /// Constructs a new test case.
        /// </summary>
        /// <param name="gridSize">
        /// <see cref="ITestCase.GridSize"/>
        /// </param>
        /// <param name="gridType">
        /// <see cref="ITestCase.GridType"/>
        /// </param>
        public TestCase(GridSizes gridSize, GridTypes gridType) {
            this.GridSize = gridSize;
            this.GridType = gridType;
        }

        /// <summary>
        /// <see cref="ITestCase.ProceedToNextShift"/>
        /// </summary>
        /// <returns>
        /// <see cref="ITestCase.ProceedToNextShift"/>
        /// </returns>
        public bool ProceedToNextShift() {
            CurrentShift = AllShifts.ElementAtOrDefault(++shiftIndex);
            return (CurrentShift != null);
        }

        /// <summary>
        /// <see cref="ITestCase.ResetShifts"/>
        /// </summary>
        public void ResetShifts() {
            shiftIndex = -1;
            CurrentShift = null;
        }

        /// <summary>
        /// <see cref="ITestCase.ScaleShifts"/>
        /// </summary>
        /// <param name="scaling">
        /// <see cref="ITestCase.ScaleShifts"/>
        /// </param>
        public void ScaleShifts(double scaling) {
            foreach (var shift in AllShifts) {
                shift.Scale(scaling);
            }
        }

        /// <summary>
        /// <see cref="ITestCase.NumberOfShifts"/>
        /// </summary>
        public int NumberOfShifts {
            get {
                return AllShifts.Count();
            }
        }

        /// <summary>
        /// <see cref="ITestCase.CurrentShift"/>
        /// </summary>
        public ShiftType CurrentShift {
            get;
            set;
        }

        /// <summary>
        /// <see cref="ITestCase.AllShifts"/>
        /// </summary>
        protected abstract IEnumerable<ShiftType> AllShifts {
            get;
        }

        /// <summary>
        /// <see cref="ITestCase.GetGrid"/>
        /// </summary>
        /// <param name="db">
        /// <see cref="ITestCase.GetGrid"/>
        /// </param>
        /// <returns>
        /// <see cref="ITestCase.GetGrid"/>
        /// </returns>
        public abstract GridCommons GetGrid(IDatabaseInfo db);

        /// <summary>
        /// <see cref="ITestCase.GridSpacing"/>
        /// </summary>
        public abstract double GridSpacing {
            get;
        }

        /// <summary>
        /// <see cref="ITestCase.Solution"/>
        /// </summary>
        public abstract double Solution {
            get;
        }

        /// <summary>
        /// <see cref="ITestCase.IntegrandDegree"/>
        /// </summary>
        public virtual int IntegrandDegree {
            get {
                return 0;
            }
        }

        /// <summary>
        /// <see cref="ITestCase.GetLevelSet"/>
        /// </summary>
        public abstract ILevelSet GetLevelSet(GridData gridData);

        /// <summary>
        /// <see cref="ITestCase.UpdateLevelSet"/>
        /// </summary>
        public abstract void UpdateLevelSet(ILevelSet levelSet);

        /// <summary>
        /// Initializes the continuous field with the value 1.
        /// </summary>
        /// <param name="input">Irrelevant</param>
        /// <param name="output">1</param>
        /// <seealso cref="ITestCase.ContinuousFieldInitialValue"/>
        public virtual void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                output[i] = 1.0;
            }
        }

        /// <summary>
        /// Initializes the discontinuous field with the value 0 in all regions
        /// associated with negative level set values
        /// </summary>
        /// <param name="input">Irrelevant</param>
        /// <param name="output">0</param>
        public virtual void JumpingFieldSpeciesAInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
        }

        /// <summary>
        /// Initializes the discontinuous field with the value 1 in all regions
        /// associated with positive level set values
        /// </summary>
        /// <param name="input">Irrelevant</param>
        /// <param name="output">1</param>
        public virtual void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output) {
            for (int i = 0; i < input.GetLength(0); i++) {
                output[i] = 1.0;
            }
        }
    }
}
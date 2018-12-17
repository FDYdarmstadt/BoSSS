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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CutCellQuadrature.TestCases {

    /// <summary>
    /// Grid solution
    /// </summary>
    public enum GridSizes {

        /// <summary>
        /// Coarsest resolution
        /// </summary>
        Tiny,

        /// <summary>
        /// Second coarsest resolution
        /// </summary>
        Small,

        /// <summary>
        /// Medium resolution
        /// </summary>
        Normal,

        /// <summary>
        /// Second finest resolution
        /// </summary>
        Large,

        /// <summary>
        /// Finest resolution
        /// </summary>
        Huge
    }

    /// <summary>
    /// The type of the grid that the test case should be performed on
    /// </summary>
    public enum GridTypes {

        /// <summary>
        /// A fully structure grid (consisting of either squares or
        /// hexahedrals)
        /// </summary>
        Structured,

        /// <summary>
        /// An unstructured triangle/tetra grid that is obtained by subdividing
        /// a structured grid)
        /// </summary>
        PseudoStructured,

        /// <summary>
        /// A fully unstructured grid provided by a grid generator (consisting
        /// of either triangles or tetras)
        /// </summary>
        Unstructured
    }

    /// <summary>
    /// A generic test case
    /// </summary>
    public interface ITestCase {

        /// <summary>
        /// Loads the next shift of the level set function-
        /// </summary>
        /// <returns>
        /// True, if there is another shift to be considered; false
        /// otherwise
        /// </returns>
        bool ProceedToNextShift();

        /// <summary>
        /// Returns to the first shift.
        /// </summary>
        void ResetShifts();

        /// <summary>
        /// Scales all shifts, see <see cref="Shift.Scale"/>
        /// </summary>
        /// <param name="scaling"></param>
        void ScaleShifts(double scaling);

        /// <summary>
        /// The total number of shifts that are defined
        /// </summary>
        int NumberOfShifts {
            get;
        }

        /// <summary>
        /// Returns the computational mesh for the current test.
        /// </summary>
        /// <param name="db">
        /// The database where the grid is stored (optionally)
        /// </param>
        /// <returns>
        /// The grid on which the test should be performed
        /// </returns>
        IGrid GetGrid(IDatabaseInfo db);

        /// <summary>
        /// The type of the grid.
        /// </summary>
        GridTypes GridType {
            get;
        }

        /// <summary>
        /// The size of the grid.
        /// </summary>
        GridSizes GridSize {
            get;
        }

        /// <summary>
        /// A measure for the size of the cells <b>on the coarsest</b> grid
        /// </summary>
        double GridSpacing {
            get;
        }

        /// <summary>
        /// The true solution of the test case
        /// </summary>
        double Solution {
            get;
        }

        /// <summary>
        /// The degree of the integrand that is considered in this test case
        /// </summary>
        int IntegrandDegree {
            get;
        }

        /// <summary>
        /// Retrieves the level set function that is associated with this test.
        /// </summary>
        /// <param name="gridData">
        /// Information about the grid.
        /// </param>
        /// <returns>
        /// The level set function associated with this test.
        /// </returns>
        ILevelSet GetLevelSet(GridData gridData);

        /// <summary>
        /// Updates the given level set according to the current shift.
        /// </summary>
        /// <param name="levelSet">
        /// The level set to be modified
        /// </param>
        void UpdateLevelSet(ILevelSet levelSet);

        /// <summary>
        /// Initial value of a continuous integrand. Required for the testing
        /// of surface integrals
        /// </summary>
        /// <param name="input">
        /// <see cref="ScalarFunction"/>
        /// </param>
        /// <param name="output">
        /// <see cref="ScalarFunction"/>
        /// </param>
        void ContinuousFieldInitialValue(MultidimensionalArray input, MultidimensionalArray output);

        /// <summary>
        /// Initial value of a discontinuous integrand in regions with negative
        /// level set values. Required for the testing of volume integrals
        /// </summary>
        /// <param name="input">
        /// <see cref="ScalarFunction"/>
        /// </param>
        /// <param name="output">
        /// <see cref="ScalarFunction"/>
        /// </param>
        void JumpingFieldSpeciesAInitialValue(MultidimensionalArray input, MultidimensionalArray output);

        /// <summary>
        /// Initial value of a discontinuous integrand in regions with positive
        /// level set values. Required for the testing of volume integrals
        /// </summary>
        /// <param name="input">
        /// <see cref="ScalarFunction"/>
        /// </param>
        /// <param name="output">
        /// <see cref="ScalarFunction"/>
        /// </param>
        void JumpingFieldSpeciesBInitialValue(MultidimensionalArray input, MultidimensionalArray output);
    }

    /// <summary>
    /// Extension methods <see cref="ITestCase"/>
    /// </summary>
    public static class ITestCaseExtensions {

        /// <summary>
        /// Retrieves the appropriate regularization polynomial.
        /// </summary>
        /// <param name="testCase">
        /// The test case under consideration.
        /// </param>
        /// <param name="numberOfVanishingMoments">
        /// The requested number of vanishing moments
        /// </param>
        /// <param name="numberOfContinuousDerivatives">
        /// The requested number of continuous derivatives.
        /// </param>
        /// <returns></returns>
        public static RegularizationPolynomoial GetPolynomial(this ITestCase testCase, int numberOfVanishingMoments, int numberOfContinuousDerivatives) {
            if (testCase is IVolumeTestCase) {
                return HeavisidePolynomial.GetPolynomial(numberOfVanishingMoments, numberOfContinuousDerivatives);
            } else if (testCase is ISurfaceTestCase) {
                return DeltaPolynomial.GetPolynomial(numberOfVanishingMoments, numberOfContinuousDerivatives);
            } else {
                throw new NotImplementedException();
            }
        }
    }

    /// <summary>
    /// Tags test cases for surface integrals
    /// </summary>
    interface ISurfaceTestCase {
    }

    /// <summary>
    /// Tags test cases for volume integrals
    /// </summary>
    interface IVolumeTestCase {
    }
}
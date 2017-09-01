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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP.Utils;
using NUnit.Framework;

namespace QuadratureAndProjectionTest {

    /// <summary>
    /// Tests for <see cref="Square"/>-related stuff.
    /// </summary>
    [TestFixture]
    public class SquareTest : QuadratueAndProjectionTest {

        /// <summary>
        /// Creates a new <see cref="Square"/>
        /// </summary>
        /// <returns></returns>
        protected override RefElement GetSimplex() {
            return Square.Instance;
        }

        /// <summary>
        /// Creates a grid containing the reference <see cref="Square"/> only
        /// </summary>
        /// <returns></returns>
        public override GridCommons GetSingleCellGrid() {
            double[] nodes = GenericBlas.Linspace(-1, 1, 2);
            return Grid2D.Cartesian2DGrid(nodes, nodes);
        }
    }
}

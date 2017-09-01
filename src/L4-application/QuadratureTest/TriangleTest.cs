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
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace QuadratureAndProjectionTest {

    /// <summary>
    /// Tests for <see cref="Triangle"/>-related stuff
    /// </summary>
    public class TriangleTest : QuadratueAndProjectionTest {

        /// <summary>
        /// Returns a new instance of <see cref="Triangle"/>
        /// </summary>
        /// <returns></returns>
        protected override RefElement GetSimplex() {
            return Triangle.Instance;
        }

        /// <summary>
        /// Returns a grid containing the reference triangle only
        /// </summary>
        /// <returns></returns>
        public override GridCommons GetSingleCellGrid() {
            MultidimensionalArray Vertices = MultidimensionalArray.Create(3, 2);
            Vertices[0, 0] = 0;
            Vertices[0, 1] = 4.0 / 3.0;
            Vertices[1, 0] = -2.0 * Math.Sqrt(3.0) / 3.0;
            Vertices[1, 1] = -2.0 / 3.0;
            Vertices[2, 0] = 2.0 * Math.Sqrt(3.0) / 3.0;
            Vertices[2, 1] = -2.0 / 3.0;

            var grd = new Grid2D(Triangle.Instance);
            grd.Cells = new Cell[1] {
                new Cell() {
                    GlobalID = 0,
                    NodeIndices = new int[] {0, 1, 2},
                    TransformationParams = Vertices,
                    Type = CellType.Triangle_3
                }
            };

            return grd;
        }
    }
}

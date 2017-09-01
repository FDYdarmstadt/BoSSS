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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using NUnit.Framework;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace QuadratureAndProjectionTest {

    /// <summary>
    /// Tests for <see cref="Tetra"/>-related components
    /// </summary>
    [TestFixture]
    public class TetraTest : QuadratueAndProjectionTest {

        /// <summary>
        /// Constructs a new <see cref="Tetra"/>
        /// </summary>
        /// <returns></returns>
        protected override RefElement GetSimplex() {
            return Tetra.Instance;
        }

        /// <summary>
        /// Creates a grid containing the reference <see cref="Tetra"/> only.
        /// </summary>
        /// <returns></returns>
        public override GridCommons GetSingleCellGrid() {
            var Simplex = Tetra.Instance;
            return new Grid3D(Simplex) {
                Cells = new Cell[] {
                    new Cell() {
                        GlobalID = 0,
                        NodeIndices = new int[] {0,1,2,3},
                        TransformationParams = Simplex.Vertices,
                        Type = CellType.Tetra_Linear
                    }
                }
            };
        }
    }
}

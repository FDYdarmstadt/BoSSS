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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP.Utils;

namespace CNS.Tests.IsentropicVortex {

    class VortexRock4 : Program<CNSControl> {

        //static void Main(string[] args) {
        //    BoSSS.Solution.Application._Main(args, false, null, () => new VortexRock4());
        //}

        protected override GridCommons CreateOrLoadGrid() {

            int NumOfCells = 4;

            double[] xnodes1 = GenericBlas.Linspace(-10, -5, NumOfCells + 1);
            double[] xnodes2 = GenericBlas.Linspace(-5, 5, 4 * NumOfCells + 1);
            double[] xnodes3 = GenericBlas.Linspace(5, 10, NumOfCells + 1);

            double[] xComplete = new double[xnodes1.Length + xnodes2.Length + xnodes3.Length - 2];
            for (int i = 0; i < xnodes1.Length; i++) {
                xComplete[i] = xnodes1[i];
            }
            for (int i = 1; i < xnodes2.Length; i++) {
                xComplete[i + xnodes1.Length - 1] = xnodes2[i];
            }
            for (int i = 1; i < xnodes3.Length; i++) {
                xComplete[i + xnodes1.Length + xnodes2.Length - 2] = xnodes3[i];
            }


            GridCommons grd = Grid2D.Cartesian2DGrid(xComplete, xComplete, CellType.Square_Linear, true, true);

            grd.Description = "2D cartesian grid 10x10 cells";

            return grd;
        }
    }
}

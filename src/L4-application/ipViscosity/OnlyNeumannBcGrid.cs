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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.ipViscosity {

    class OnlyNeumannBcGrid : ITestGrid {

        public GridCommons GetGrid() {

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 7), GenericBlas.Linspace(-2, 2, 5));

            grd.EdgeTagNames.Add(1, "Pressure_Outlet");

            grd.DefineEdgeTags(delegate(double[] _X) {

                return 1;
            });


            return grd;
        }


        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            // values are overriden anyway
            return new Dictionary<string, AppControl.BoundaryValueCollection>() {
                { "Pressure_Outlet", new AppControl.BoundaryValueCollection() },
            };
        }
    }
}

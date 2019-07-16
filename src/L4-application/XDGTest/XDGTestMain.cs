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
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using ilPSP.Utils;
using ilPSP.Tracing;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Application.XDGTest {

    /// <summary>
    /// Some basic tests for the XDG framework, e.g. DG projection and extrapolation under level-set movement
    /// </summary>
    class XDGTestMain : BoSSS.Solution.Application {
        static void Main(string[] args) {
            
            BoSSS.Solution.Application._Main(args, true, delegate() {
                var p = new XDGTestMain();
                //p.passiveIo = true;
                return p;
            });
        }

        /// <summary>
        /// pressure computed in projection step
        /// </summary>
        XDGField Pressure;

        /// <summary>
        /// the level set
        /// </summary>
        LevelSet LevSet;

        protected override void CreateFields() {
            this.LevSet = new LevelSet(new Basis(this.gridData, 2), "Phi");
            this.LsTrk = new LevelSetTracker((GridData) this.gridData, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, LevSet);
            Pressure = new XDGField(new XDGBasis(this.LsTrk, 2), "Pressure");
        }

        protected override IGrid CreateOrLoadGrid() {
            var xNodes = GenericBlas.Linspace(-0.33333, 0.666667, 7);
            //var yNodes = GenericBlas.Linspace(-1, 1, 2);
            //var xNodes = GenericBlas.Linspace(-2, 2, 25);
            var yNodes = GenericBlas.Linspace(-1, 1, 13);
            var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
            return grd;
        }

        protected override void SetInitial() {
            //this.LevSet.ProjectField((_2D)((x, y) => (x - 0.1)));
            this.LevSet.ProjectField((_2D)((x, y) => ((x - 0.83) / 0.8).Pow2() + (y / 0.8).Pow2() - 1.0));
            this.LsTrk.UpdateTracker();
            this.LsTrk.PushStacks();
            this.Pressure.ProjectField((_2D)((x, y) => 1 - y * y));

            //PlotCurrentState(0, 0, 3);
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
        }

        internal double AutoExtrapolationErr = double.MinValue;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            dt = 1.0;
            base.NoOfTimesteps = 1;

            Console.WriteLine("Timestep #{0}, dt = {1} ...", TimestepNo, dt);

            // test the auto-extrapolation
            // ---------------------------

            double x0 = 0.17 * (phystime + dt);
            this.Pressure.UpdateBehaviour = BehaveUnder_LevSetMoovement.AutoExtrapolate;
            this.LevSet.ProjectField((_2D)((x, y) => ((x - 0.83 - x0) / 0.8).Pow2() + (y / 0.8).Pow2() - 1.0));
            this.LsTrk.UpdateTracker();


            var RefPressure = new XDGField(this.Pressure.Basis);
            RefPressure.ProjectField((_2D)((x, y) => 1 - y * y));
            RefPressure.Acc(-1.0, Pressure);

            AutoExtrapolationErr = RefPressure.L2Norm();
            Console.WriteLine("Error of extrapolation: " + AutoExtrapolationErr);



           // PlotCurrentState(phystime + dt, TimestepNo);


            Console.WriteLine("done.");

            return dt;
        }


        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { this.Pressure, this.LevSet }, "XNSE_prj" + timestepNo, physTime, superSampling);
        }


    }
}

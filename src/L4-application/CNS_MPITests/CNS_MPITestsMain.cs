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

using BoSSS.Foundation.IO;
using CNS;
using System;

namespace BoSSS.Solution.Template {

    /// <summary>
    /// Main application class.
    /// </summary>
    public class TemplateMain : Program<CNSControl> {

        ///// <summary>
        ///// Application entry point.
        ///// </summary>
        //static void Main(string[] args) {
        //    Application<CNSControl>._Main(args, true, null, delegate () {
        //        return new TemplateMain();
        //    });
        //}

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            throw new NotImplementedException();
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase loadBalancingData) {
            base.CreateEquationsAndSolvers(loadBalancingData);
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            throw new NotImplementedException();
        }

    }
}

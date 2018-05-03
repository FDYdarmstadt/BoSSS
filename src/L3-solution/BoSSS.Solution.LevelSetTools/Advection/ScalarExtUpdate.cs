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
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using ilPSP;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.LevelSetTools.Advection {
    /// <summary>
    /// A class for obtaining the time integrated spatial discretization of the level set advection equation
    /// </summary>
    public class ScalarExtUpdate : ILevelSetAdvection {

        Extender VelocityExtender;
        SinglePhaseField Extension;
        SinglePhaseField LevelSet;
        bool nearfield;

        public ScalarExtUpdate(LevelSetTracker LSTrk, SinglePhaseField LevelSet, VectorField<SinglePhaseField> LevelSetGradient, VectorField<SinglePhaseField> Velocity, EllipticExtVelAlgoControl Control, out SinglePhaseField Extension, bool nearfield) {

            this.LevelSet = LevelSet;
            this.nearfield = nearfield;

            int D = LSTrk.GridDat.SpatialDimension;
            double PenaltyBase = Control.PenaltyMultiplierInterface * ((double)((LevelSet.Basis.Degree + 1) * (LevelSet.Basis.Degree + D))) / ((double)D);

            ILevelSetForm InterfaceFlux = new ScalarVelocityInterfaceForm(PenaltyBase, LSTrk);

            List<DGField> Paramlist = new List<DGField> { };
            Paramlist.AddRange(Velocity);

            Extension = new SinglePhaseField(LevelSet.Basis,"ExtensionVelocity");
            this.Extension = Extension;
            VelocityExtender = new Extender(Extension, LSTrk, InterfaceFlux, Paramlist, LevelSetGradient, Control);
            VelocityExtender.ConstructExtension(nearfield:nearfield);
        }

        public void Advect(double dt) {
            VelocityExtender.ConstructExtension(nearfield:nearfield);
            LevelSet.Acc(-dt, Extension);

        }

        public void FinishTimeStep() {
            // Do nothing so far
        }
    }
}
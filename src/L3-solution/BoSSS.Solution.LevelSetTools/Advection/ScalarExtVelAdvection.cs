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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.EllipticExtension;

namespace BoSSS.Solution.LevelSetTools.Advection {
    public class ScalarExtVelAdvection : ILevelSetAdvection {

        Extender VelocityExtender;
        ScalarVelocityAdvection Advection;
        SinglePhaseField Extension;

        bool nearfield;

        public ScalarExtVelAdvection(LevelSetTracker LSTrk, SinglePhaseField LevelSet, VectorField<SinglePhaseField> LevelSetGradient, VectorField<SinglePhaseField> Velocity, EllipticExtVelAlgoControl Control, bool nearfield, out SinglePhaseField Extension) {

            int D = LSTrk.GridDat.SpatialDimension;
            double PenaltyBase = Control.PenaltyMultiplierInterface * ((double)((LevelSet.Basis.Degree + 1) * (LevelSet.Basis.Degree + D))) / ((double)D);

            ILevelSetForm InterfaceFlux = new ScalarVelocityInterfaceForm(PenaltyBase, LSTrk);

            this.nearfield = nearfield;

            List<DGField> Paramlist = new List<DGField> { };
            Paramlist.AddRange(Velocity);

            Extension = new SinglePhaseField(Velocity[0].Basis, "ExtensionVelocity");

            VelocityExtender = new Extender(Extension, LSTrk, InterfaceFlux, Paramlist, LevelSetGradient, Control);
            Advection = new ScalarVelocityAdvection(LSTrk, LevelSet, Extension, null, nearfield);
        }



        public void Advect(double dt) {
            VelocityExtender.ConstructExtension(nearfield:nearfield);
            Advection.Advect(dt);
        }

        public void FinishTimeStep() {
            // Do nothing so far
        }
    }
}

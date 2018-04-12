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
    public class VectorExtVelAdvection : ILevelSetAdvection {

        Extender[] VelocityExtender;
        VectorVelocityLevelSetAdvection Advection;
        VectorField<SinglePhaseField> Velocity;
        SinglePhaseField LevelSet;
        bool nearfield;
        VectorFieldHistory<SinglePhaseField> VectorExtensionHistory;
        VectorField<SinglePhaseField> VectorExtension;
        VectorField<SinglePhaseField> LevelSetGradient;

        int D;
        bool UseImplicitTimeStepping;
        LevelSetTracker LSTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="LSTrk"></param>
        /// <param name="LevelSet"></param>
        /// <param name="LevelSetGradient"></param>
        /// <param name="Velocity"></param>
        /// <param name="Control"></param>
        /// <param name="bcMap"></param>
        /// <param name="BDForder">Quick Hack: -1 = CrankNicholson, 0= Explicit, >0: BDF</param>
        /// <param name="VectorExtension"></param>
        /// <param name="nearfield"></param>
        public VectorExtVelAdvection(LevelSetTracker LSTrk, SinglePhaseField LevelSet,VectorField<SinglePhaseField> LevelSetGradient, VectorField<SinglePhaseField> Velocity, EllipticExtVelAlgoControl Control, NSECommon.IncompressibleBoundaryCondMap bcMap, int BDForder, out VectorField<SinglePhaseField> VectorExtension, SubGrid subGrid = null) {

            D = LSTrk.GridDat.SpatialDimension;
            this.LevelSetGradient = LevelSetGradient;
            this.UseImplicitTimeStepping = (BDForder != 0);
            this.nearfield = subGrid !=null ;
            this.Velocity = Velocity;
            this.LSTrk = LSTrk;
            this.LevelSet = LevelSet;

            

            VectorExtensionHistory = new VectorFieldHistory<SinglePhaseField>(new VectorField<SinglePhaseField>(D, Velocity[0].Basis, "ExtVel", SinglePhaseField.Factory));
            VectorExtension = VectorExtensionHistory.Current;
            this.VectorExtension = VectorExtension;
            


            double PenaltyBase = Control.PenaltyMultiplierInterface * ((double)((LevelSet.Basis.Degree + 1) * (LevelSet.Basis.Degree + D))) / ((double)D);
            ILevelSetComponent InterfaceFlux = new SingleComponentInterfaceForm(PenaltyBase, LSTrk);


            VelocityExtender = new Extender[D];
            for (int d = 0; d < D; d++) {
                VelocityExtender[d] = new Extender( VectorExtension[d] , LSTrk, InterfaceFlux, new List<DGField> {Velocity[d] }, LevelSetGradient, Control);
                VelocityExtender[d].ConstructExtension(new List <DGField> { Velocity[d] }, Control.subGridRestriction);
            }


            if (!this.UseImplicitTimeStepping) {
                VectorExtensionHistory.IncreaseHistoryLength(1);
                VectorExtensionHistory.Push();

                Advection = new ExplicitNonconservativeAdvection(LevelSet, VectorExtensionHistory, bcMap);
            }
            else {
                    Advection = new BDFNonconservativeAdvection(LevelSet, VectorExtension, bcMap, BDForder, subGrid, false);
            }
         }


        
        public void Advect(double dt) {
            VectorExtension.Clear();
            for (int d = 0; d < D; d++) {
                VelocityExtender[d].ConstructExtension(new List <DGField> { Velocity[d] }, this.nearfield);
            }

            if (!this.UseImplicitTimeStepping) VectorExtensionHistory.Push();

            Advection.Advect(dt);
        }

        public void FinishTimeStep() {
            Advection.FinishTimeStep();
        }
    }
}

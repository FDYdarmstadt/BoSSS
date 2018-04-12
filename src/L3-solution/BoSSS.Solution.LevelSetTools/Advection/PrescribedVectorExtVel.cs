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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.LevelSetTools.Advection
{
    public class PrescribedVectorExtVel : ILevelSetAdvection
    {
        SinglePhaseField LevelSet;
        VectorField<SinglePhaseField> Velocity;
        VectorFieldHistory<SinglePhaseField> VelocityHistory;
        double t = 0;
        bool nearfield;
        Func<double[], double, double[]> AnalyticalExtVel;
        ExplicitNonconservativeAdvection explicitNonconservativeAdvection;




        public PrescribedVectorExtVel(SinglePhaseField LevelSet, VectorField<SinglePhaseField> Velocity, bool nearfield, NSECommon.IncompressibleBoundaryCondMap bcMap, Func<double[], double, double[]> AnalyticalExtVel) {
            this.nearfield = nearfield;
            this.LevelSet = LevelSet;
            this.Velocity = Velocity;
            this.AnalyticalExtVel = AnalyticalExtVel;

            VelocityHistory = new VectorFieldHistory<SinglePhaseField>(Velocity);
            VelocityHistory.IncreaseHistoryLength(1);
            VelocityHistory.Push();
            explicitNonconservativeAdvection = new ExplicitNonconservativeAdvection(LevelSet, VelocityHistory, bcMap, false);
        }
        
        public void Advect(double dt){
            Func<double[], double> func;
            for (int d=0; d<Velocity.Dim; d++) {
                func = new Func<double[], double>(X => AnalyticalExtVel(X, t)[d]);
                VelocityHistory.Current[d].Clear();
                VelocityHistory.Current[d].ProjectField(func);
            }
            explicitNonconservativeAdvection.Advect(dt);
            VelocityHistory.Push();
            t += dt;

        }

        public void FinishTimeStep()
        {
            throw new NotImplementedException();
        }
    }
}

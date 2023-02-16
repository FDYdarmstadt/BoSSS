using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class GhostPenalty : BulkEquation {

        string spcName;

        string codomain; 

        public GhostPenalty(string spcName, int d, int D, double scale) {
            this.spcName = spcName;
            this.codomain = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);

            for(int i = 0; i < D; ++i) {
                string velocity = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[i];
                AddVariableNames(velocity);
                var divergence1 = new EdgeGradientPenaltyForm(spcName, velocity, scale);
                AddComponent(divergence1);
            }
        }

        public override string SpeciesName => spcName;

        public override double MassScale => 0.0;

        public override string CodomainName => codomain;
    }
}

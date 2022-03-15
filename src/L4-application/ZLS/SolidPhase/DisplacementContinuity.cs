using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class DisplacementContinuity : BulkEquation {

        internal static bool ContinuityStabilization = false;

        string spcName;

        public DisplacementContinuity(string spcName, int D, Solid Material) {
            this.spcName = spcName;
            for(int i = 0; i < D; ++i) {
                string displacement = ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[i];
                AddVariableNames(displacement);
                var divergence = new Divergence(spcName, displacement, i, 1.0);
                AddComponent(divergence);
            }
        }

        public override string SpeciesName => spcName;

        public override double MassScale => 0.0;

        public override string CodomainName => ZwoLevelSetSolver.EquationNames.DisplacementContinuity;
    }
}

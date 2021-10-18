using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class PressurePenalty : BulkEquation {

        string speciesName;

        //Solid material;

        string codomainName;

        public PressurePenalty(string speciesName, double scale) {
            this.speciesName = speciesName;
            //this.material = material;
            this.codomainName = BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);

            var pressurePenalty = new EdgePenaltyForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.Pressure, scale);
            AddComponent(pressurePenalty);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 0;

        public override string CodomainName => codomainName;

    }
}

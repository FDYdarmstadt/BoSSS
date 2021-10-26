using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class Dummy : BulkEquation {

        string speciesName;
        string codomainName;

        public Dummy(string speciesName, string domainName, string codomainName) {
            this.speciesName = speciesName;
            this.codomainName = codomainName;
            AddVariableNames(domainName);
            var variable = new MultiPhaseVariableSource(speciesName, domainName);
            AddComponent(variable);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 0.0;

        public override string CodomainName => codomainName;
    }
}

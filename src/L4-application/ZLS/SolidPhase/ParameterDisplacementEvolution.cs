using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon.Operator;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.SolidPhase {
    class ParameterDisplacementEvolution : BulkEquation {

        string speciesName;

        string codomainName;

        public ParameterDisplacementEvolution(string speciesName, int d, int D, double artificialViscosity) {
            this.speciesName = speciesName;
            this.codomainName = EquationNames.DisplacementEvolutionComponent(d);
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));

            var convection = new LinearTransportForm(speciesName, ZwoLevelSetSolver.VariableNames.Displacement0Vector(D), 
                ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, D, 1.0);
            AddComponent(convection);
            AddParameter(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D));

            var source = new MultiPhaseSource(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D)[d], speciesName, -1.0);
            AddComponent(source);

            if(artificialViscosity != 0) {
                // we should not add the SIP form if it is not intended at all, i.e. if 'artificialViscosity == 0';
                // since evaluation of SIP forms is quite costly; 
                AddComponent(new SIPForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, artificialViscosity));
            }
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }
}

using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ContactLine {
    class EquilibriumContactLine : SurfaceEquation {
        string codomainName;

        public override string FirstSpeciesName => "A";

        public override string SecondSpeciesName => "B";

        public override string CodomainName => codomainName;

        public EquilibriumContactLine(int d, int D, double beta, double theta) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddContactLineComponent(new LaplaceBeltramiEquilibriumForm(d, D));

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)
                .Cat(BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(ZwoLevelSetSolver.VariableNames.SolidLevelSetCG, BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)))
                .Cat(BoSSS.Solution.NSECommon.VariableNames.MaxSigma));
        }
    }

    class SlipContactLine : SurfaceEquation {
        string codomainName;

        public override string FirstSpeciesName => "A";

        public override string SecondSpeciesName => "B";

        public override string CodomainName => codomainName;

        public SlipContactLine(int d, int D, double beta, double theta) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddContactLineComponent(new LaplaceBeltramiSlipForm(d, D));

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)
                .Cat(BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(ZwoLevelSetSolver.VariableNames.SolidLevelSetCG, BoSSS.Solution.NSECommon.VariableNames.NormalVector(D)))
                .Cat(BoSSS.Solution.NSECommon.VariableNames.MaxSigma));
        }
    }
}

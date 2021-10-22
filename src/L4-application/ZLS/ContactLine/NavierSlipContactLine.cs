using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using NSEVariableNames = BoSSS.Solution.NSECommon.VariableNames;
using ZLSVariableNames = ZwoLevelSetSolver.VariableNames;

namespace ZwoLevelSetSolver.ContactLine {

    class NavierSlipContactLine : SurfaceEquation {

        string codomainName;

        public override string FirstSpeciesName => "A";

        public override string SecondSpeciesName => "B";

        public override string CodomainName => codomainName;

        public NavierSlipContactLine(int d, int D, double beta, double theta) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddContactLineComponent(new NavierSlipLinearContactLineForm(d, D, beta, theta));
            AddContactLineComponent(new FreeSlipContactLineForm(d, D));

            AddVariableNames(NSEVariableNames.VelocityVector(D));
            AddParameter(NSEVariableNames.NormalVector(D)
                .Cat(NSEVariableNames.AsLevelSetVariable( ZLSVariableNames.SolidLevelSetCG, NSEVariableNames.NormalVector(D)))
                .Cat(NSEVariableNames.MaxSigma));
        }
    }
}

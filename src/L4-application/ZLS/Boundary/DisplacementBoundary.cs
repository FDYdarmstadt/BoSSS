using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class DisplacementBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public DisplacementBoundary(LevelSetTracker LsTrkr, string fluidSpecies, string solidSpecies, int d, int D, double artificialViscosity) {
            codomainName = ZwoLevelSetSolver.EquationNames.DisplacementEvolutionComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddComponent(new NonLinearSolidConvectionForm(ZwoLevelSetSolver.VariableNames.DisplacementVector(D), BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), 1.0, d, 1, fluidSpecies, solidSpecies));
            AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, D, 1, artificialViscosity));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

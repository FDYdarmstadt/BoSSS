using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class LinearDisplacementBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public LinearDisplacementBoundary(LevelSetTracker LsTrkr, string fluidSpecies, string solidSpecies, int d, int D, double artificialViscosity) {
            codomainName = ZwoLevelSetSolver.EquationNames.DisplacementEvolutionComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddComponent(new DisplacementBoundaryConvectionForm(ZwoLevelSetSolver.VariableNames.DisplacementVector(D), 1.0,d, D, 1, fluidSpecies, solidSpecies));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));

            AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, D, 1, artificialViscosity));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

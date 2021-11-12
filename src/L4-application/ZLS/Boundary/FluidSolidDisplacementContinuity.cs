using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class FluidSolidDisplacementContinuity : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public FluidSolidDisplacementContinuity(string fluidSpecies, string solidSpecies, int D) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality

            for(int i = 0; i < D; ++i) {
                string velocity = ZwoLevelSetSolver.VariableNames.DisplacementComponent(i);
                AddVariableNames(velocity);
                AddComponent(new BoundaryDivergenceDisplacementForm(i, 1, fluidSpecies, solidSpecies));
            }
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }

}

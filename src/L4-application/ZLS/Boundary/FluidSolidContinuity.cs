using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {
    class FluidSolidContinuity : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public FluidSolidContinuity(string fluidSpecies, string solidSpecies, int D) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality
            
            for(int i = 0; i < D; ++i) {
                string velocity = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[i];
                AddVariableNames(velocity);
                AddComponent(new BoundaryDivergenceForm( velocity, i, 1, fluidSpecies, solidSpecies));
            }
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

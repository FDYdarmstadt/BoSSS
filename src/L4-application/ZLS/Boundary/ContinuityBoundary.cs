using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.Boundary {
    class ContinuityBoundary : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public ContinuityBoundary(string fluidSpecies, string solidSpecies, int D) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality

            for(int i = 0; i < D; ++i) {
                string velocity = ZwoLevelSetSolver.VariableNames.DisplacementComponent(i);
                //string velocity = BoSSS.Solution.NSECommon.VariableNames.Velocity_d(i);
                AddVariableNames(velocity);
                AddComponent(new BoundaryDivergenceForm(i, 1, fluidSpecies, solidSpecies, 1));
            }
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }

}

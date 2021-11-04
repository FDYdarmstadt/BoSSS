using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.ContactLine;

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
                
                //string variableName = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[i];
                //string variableName = ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[i];
                //AddSurfaceComponent(new SurfacePenaltyForm(variableName, -100));
            }
            //AddSurfaceComponent(new SurfacePenaltyForm(BoSSS.Solution.NSECommon.VariableNames.Pressure, 100));
            //AddSurfaceComponent(new SurfaceGradientPenaltyForm(BoSSS.Solution.NSECommon.VariableNames.Pressure, D, 10));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

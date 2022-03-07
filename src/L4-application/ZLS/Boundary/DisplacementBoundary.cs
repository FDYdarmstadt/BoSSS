using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.ContactLine;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.Boundary {
    
    /// <summary>
    /// Solid Fluid Interface for Displacement Transport Equation
    /// </summary>
    class DisplacementBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;


        public DisplacementBoundary(LevelSetTracker LsTrkr, string fluidSpecies, string solidSpecies, int d, int D, double artificialViscosity, double fluidViscosity, Solid material) {
            codomainName = ZwoLevelSetSolver.EquationNames.DisplacementEvolutionComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddComponent(new InnerNonLinearSolidConvectionForm(ZwoLevelSetSolver.VariableNames.DisplacementVector(D), 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), 
                1.0, d, 1, 
                fluidSpecies, 
                solidSpecies));
            if (artificialViscosity != 0.0)
            {
                AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, d, D, 1, artificialViscosity, fluidViscosity, material.Viscosity));
            }
            //AddSurfaceComponent(new SurfacePenaltyForm(ZwoLevelSetSolver.VariableNames.DisplacementComponent(d), -100));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }

    class ExtensionDisplacementBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public ExtensionDisplacementBoundary(LevelSetTracker LsTrkr, string fluidSpecies, string solidSpecies, int d, int D, double extensionViscosity, double artificialViscosity) {
            codomainName = ZwoLevelSetSolver.EquationNames.DisplacementEvolutionComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddComponent(new NonLinearSolidConvectionForm(ZwoLevelSetSolver.VariableNames.DisplacementVector(D),
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D),
                1.0, d, 1,
                fluidSpecies,
                solidSpecies));
            if(artificialViscosity != 0.0) {
                AddComponent(new ExtensionSolidTensionForm(
                    fluidSpecies, solidSpecies, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, D, 1, extensionViscosity, artificialViscosity));
            }
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

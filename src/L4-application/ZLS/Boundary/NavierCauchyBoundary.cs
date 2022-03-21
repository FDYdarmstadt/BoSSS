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
    class NavierCauchyBoundary : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public NavierCauchyBoundary(string fluidSpecies, string solidSpecies, int d, int D, Solid material, double rho_fluid, double viscosity) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));

            //Stress equality
            //AddComponent(new SolidLinearIncompressibleNeoHookeanBoundaryForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            AddComponent(new NeoHookeanNeumannForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            //AddComponent(new NonLinearNeoHookeanNeumannForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            //AddComponent(new SlipSolidLinearIncompressibleNeoHookeanBoundaryForm(fluidSpecies, solidSpecies, d, 1, viscosity,material.Viscosity, material.Lame2));

            
            //Penalty coupling
            AddComponent(new NoSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Viscosity));
            //AddComponent(new NavierSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2, 0.1));

            /*
            AddComponent(new NonLinearSolidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), material.Density, D, 1, fluidSpecies, solidSpecies));
            
            AddComponent(new NonLinearFluidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), rho_fluid, D, 1, fluidSpecies, solidSpecies));
            //*/
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }

    class ExtensionNavierCauchyBoundary : NavierCauchyBoundary {

        public ExtensionNavierCauchyBoundary(string fluidSpecies, string solidSpecies, int d, int D, Solid material, double rho_fluid, double viscosity) 
            : base(fluidSpecies, solidSpecies, d, D, material, rho_fluid, viscosity){
            AddComponent(new DisplacementPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, material.Lame2));
        }
    }
}

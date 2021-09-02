using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
            //AddComponent(new IncompressibleNeoHookeanBoundaryForm(fluidSpecies, solidSpecies, 1, d, viscosity, material.Lame2));
            AddComponent(new SolidLinearIncompressibleNeoHookeanBoundaryForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Lame2));

            //Penalty coupling
            AddComponent(new NoSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2));

            AddComponent(new NonLinearSolidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), material.Density, D, 1, fluidSpecies, solidSpecies));
            AddComponent(new NonLinearFluidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), rho_fluid, D, 1, fluidSpecies, solidSpecies));
            
            AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, D, 1, material.Viscosity));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

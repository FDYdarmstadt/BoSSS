using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using HFSISolver.ContactLine;
using HFSISolver.SolidPhase;

namespace HFSISolver.Boundary {
    class NavierCauchyBoundary : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public NavierCauchyBoundary(string fluidSpecies, string solidSpecies, int d, int D, 
            Solid material, double rho_fluid, double viscosity) {

            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(HFSISolver.VariableNames.DisplacementVector(D));

            //Stress equality
            //AddComponent(new NeoHookeanNeumannForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            //AddComponent(new NonLinearNeoHookeanNeumannForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            AddComponent(new NeoHookeanBoundaryForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Viscosity, material.Lame2));

            //Transport terms 
            AddComponent(new NonLinearSolidConvectionForm(
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D),
                material.Density, rho_fluid, d, 1, fluidSpecies, solidSpecies));

            //Penalty coupling
            double maxViscosity = Math.Max(viscosity, material.Viscosity);
            AddComponent(new NoSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, maxViscosity, maxViscosity));
            //AddComponent(new NavierSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2, 0.1));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

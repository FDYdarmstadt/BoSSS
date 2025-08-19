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

        public NavierCauchyBoundary(string fluidSpecies, string solidSpecies, int d, int D, 
            Solid material, double rho_fluid, double viscosity, double sliplength) {

            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));

            //Stress equality
            //AddComponent(new NeoHookeanNeumannForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            //AddComponent(new NonLinearNeoHookeanNeumannForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Viscosity, material.Lame2));
            AddComponent(new NeoHookeanBoundaryForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Viscosity, material.Lame2));

            //Transport therms 
            AddComponent(new NonLinearSolidConvectionForm(
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D),
                material.Density, rho_fluid, d, 1, fluidSpecies, solidSpecies));

            //Penalty coupling
            double maxViscosity = Math.Max(viscosity, material.Viscosity);
            if(sliplength == 0.0) {
                Console.Write("using no slip. ");
                AddComponent(new NoSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, maxViscosity, maxViscosity));
            } else if(double.IsInfinity(sliplength)) {
                Console.Write("using slip. ");
                AddComponent(new SlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, maxViscosity, maxViscosity));
            } else {
                Console.Write("using navier slip. ");
                AddComponent(new NavierSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2, sliplength));
            }
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

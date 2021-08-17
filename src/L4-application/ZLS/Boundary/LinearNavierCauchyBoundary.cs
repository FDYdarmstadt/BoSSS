using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.Boundary {
    class LinearNavierCauchyBoundary : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public LinearNavierCauchyBoundary(string fluidSpecies, string solidSpecies, int d, int D, Solid material, double rho_fluid, double viscosity, double artificialViscosity) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));

            AddComponent(new SolidLinearIncompressibleNeoHookeanBoundaryForm(fluidSpecies, solidSpecies, d, 1, viscosity, material.Lame2));

            //Penalty coupling
            AddComponent(new NoSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2));

            AddComponent(new SolidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], material.Density, D, 1, fluidSpecies, solidSpecies));
            AddComponent(new FluidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], rho_fluid, D, 1, fluidSpecies, solidSpecies));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));

            AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, D, 1, artificialViscosity));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}


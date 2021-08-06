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


            if(d == 0) {
                AddComponent(new SolidLinearIncompressibleNeoHookeanBoundaryFormX(fluidSpecies, solidSpecies, 1, viscosity, material.Lame2));
                AddComponent(new FluidLinearIncompressibleNeoHookeanBoundaryFormX(fluidSpecies, solidSpecies, 1, viscosity, material.Lame2));
            } else if(d == 1) {
                AddComponent(new SolidLinearIncompressibleNeoHookeanBoundaryFormY(fluidSpecies, solidSpecies, 1, viscosity, material.Lame2));
                AddComponent(new FluidLinearIncompressibleNeoHookeanBoundaryFormY(fluidSpecies, solidSpecies, 1, viscosity, material.Lame2));
            } else {
                throw new Exception("Spatial Dimension not supported.");
            }

            //Penalty coupling
            AddComponent(new NavierSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2, 0.1));

            AddComponent(new SolidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], material.Density, D, 1, fluidSpecies, solidSpecies));
            AddComponent(new FluidMomentumConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], rho_fluid, D, 1, fluidSpecies, solidSpecies));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));


            AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, d, D, 1, artificialViscosity));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}


using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.Boundary {
    class SimplifiedMomentumBoundary : SurfaceEquation {

        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public SimplifiedMomentumBoundary(string fluidSpecies, string solidSpecies, int d, int D, Solid material, double viscosity) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            //Stress equality
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));

            //Consistency
            //AddComponent(new SolidBoundarySIPForm(fluidID, solidID, d, D, 1, material.Lame2));
            //AddComponent(new SolidBoundaryPressureForm(d, 1, fluidID, solidID));
            AddComponent(new SolidTensionForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2));

            //AddComponent(new FluidBoundaryViscosityForm(fluidID, solidID, d, D, 1, viscosity));
            //AddComponent(new FluidBoundaryPressureForm(d, D, 1, fluidID, solidID));
            AddComponent(new FluidTensionForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2));

            //Penalty coupling
            AddComponent(new NoSlipVelocityPenaltyForm(fluidSpecies, solidSpecies, d, D, 1, viscosity, material.Lame2));
            //AddComponent(new DisplacementPenaltyForm(fluidID, solidID, d, D, 1));
            //AddComponent(new SolidTensionPenaltyForm(fluidID, solidID, d, D, 1, viscosity, material.Lame2));
            //AddComponent(new FluidTensionPenaltyForm(fluidID, solidID, d, D, 1, viscosity, material.Lame2));

            //AddComponent(new ConvectionForm(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], material.Density, D, 1, fluidID, solidID));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

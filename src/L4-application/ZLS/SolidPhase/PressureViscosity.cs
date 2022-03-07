using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class PressureViscosity : BulkEquation {
        string speciesName;

        //Solid material;

        string codomainName;

        public PressureViscosity(string speciesName, double scale) {
            this.speciesName = speciesName;
            //this.material = material;
            this.codomainName = BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);

            var pressureViscosity = new OpenSIPForm(SpeciesName, new string[] { BoSSS.Solution.NSECommon.VariableNames.Pressure }, 0, scale);
            AddComponent(pressureViscosity);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 0;

        public override string CodomainName => codomainName;
    }

    class PressureViscosityBoundary : SurfaceEquation {
        string fluidSpecies;
        string solidSpecies;
        string codomainName;

        public PressureViscosityBoundary(string fluidSpecies, string solidSpecies, int D, double viscosity) {
            codomainName = BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;

            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddComponent(new SIPBoundaryForm(fluidSpecies, solidSpecies, BoSSS.Solution.NSECommon.VariableNames.Pressure, D, viscosity, 1));
        }

        public override string FirstSpeciesName => fluidSpecies;

        public override string SecondSpeciesName => solidSpecies;

        public override string CodomainName => codomainName;
    }
}

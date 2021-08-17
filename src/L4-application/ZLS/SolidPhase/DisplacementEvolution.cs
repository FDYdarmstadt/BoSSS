using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon.Operator;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.SolidPhase {
    class DisplacementEvolution : BulkEquation {

        string speciesName;

        string codomainName; 

        public DisplacementEvolution(string speciesName, int d, int D, double artificialViscosity) {
            this.speciesName = speciesName;
            this.codomainName = EquationNames.DisplacementEvolutionComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            
            var convection = new NonLinearConvectionForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[d], BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), D, 1.0);
            AddComponent(convection);

            AddComponent(new SIPForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, artificialViscosity));

            var source = new MultiPhaseVariableSource(speciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], -1.0);
            AddComponent(source);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }

    class ParameterDisplacementEvolution : BulkEquation {

        string speciesName;

        string codomainName;

        public ParameterDisplacementEvolution(string speciesName, int d, int D) {
            this.speciesName = speciesName;
            this.codomainName = EquationNames.DisplacementEvolutionComponent(d);
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));


            var convection = new LinearConvectionForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[d], D, 1.0);
            AddComponent(convection);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D));
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));

            var source = new MultiPhaseSource(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d], speciesName, -1.0);
            AddComponent(source);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }
}

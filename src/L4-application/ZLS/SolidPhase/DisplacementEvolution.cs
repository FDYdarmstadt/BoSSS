using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon;
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

        public DisplacementEvolution(string speciesName, int d, int D, double artificialViscosity, IncompressibleMultiphaseBoundaryCondMap boundaryMap) {
            this.speciesName = speciesName;
            this.codomainName = EquationNames.DisplacementEvolutionComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));

            
            var convection = new NonLinearConvectionForm(speciesName, 
                ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[d], 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), 
                d, 1.0, boundaryMap);
            AddComponent(convection);

            if(artificialViscosity != 0) {
                // we should not add the SIP form if it is not intended at all, i.e. if 'artificialViscosity == 0';
                // since evaluation of SIP forms is quite costly; 
                AddComponent(new SIPForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, artificialViscosity, boundaryMap));
            }

            var source = new MultiPhaseVariableSource(speciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], -1.0);
            //var source = new MultiPhaseSource(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D)[d], SpeciesName, -1.0);
            //AddParameter(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D)[d]);
            AddComponent(source);

            //AddComponent(new EdgePenaltyForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[d], 1));
            //AddComponent(new EdgeGradientPenaltyForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[d], -1));
            //Console.WriteLine("Displacement evo deakt");
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

            var convection = new LinearTransportForm(speciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, D, 1.0);
            AddComponent(convection);
            AddParameter(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D));

            var source = new MultiPhaseSource(ZwoLevelSetSolver.VariableNames.Displacement0Vector(D)[d], speciesName, -1.0);
            AddComponent(source);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }
}

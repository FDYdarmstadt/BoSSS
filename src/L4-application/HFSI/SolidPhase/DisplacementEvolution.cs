using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using HFSISolver.SolidPhase;

namespace HFSISolver.SolidPhase {
    class DisplacementEvolution : BulkEquation {

        string speciesName;

        string codomainName;

        public DisplacementEvolution(string speciesName, Solid material, int d, int D, IncompressibleMultiphaseBoundaryCondMap boundaryMap) {
            this.speciesName = speciesName;
            this.codomainName = EquationNames.DisplacementEvolutionComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(HFSISolver.VariableNames.DisplacementVector(D));
            var convection = new SourceConvectionForm(speciesName,
                HFSISolver.VariableNames.DisplacementVector(D)[d],
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D),
                1.0);
            AddComponent(convection);
            var source = new MultiPhaseVariableSource(speciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], -1.0);
            AddComponent(source);
            //AddComponent(new EdgePenaltyForm1(SpeciesName, HFSISolver.VariableNames.DisplacementVector(D)[d], material.Lame2 ));
            //AddComponent(new BoundaryEdgePenaltyForm(SpeciesName, HFSISolver.VariableNames.DisplacementVector(D)[d], material.Lame2));
            //AddComponent(new EdgePenaltyForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], - material.Viscosity));
            //AddComponent(new BoundaryEdgePenaltyForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], material.Viscosity));
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;

        public override string CodomainName => codomainName;
    }
}
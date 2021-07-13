using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class SimplifiedNavierCauchy : BulkEquation {

        string speciesName;

        Solid material;

        string codomainName;

        public SimplifiedNavierCauchy(string speciesName, Solid material, int d, int D) {
            this.speciesName = speciesName;
            this.material = material;
            this.codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);

            var convection = new LinearConvectionForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], D, material.Density);
            AddComponent(convection);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d]);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]);
            
            var elasticTension = new SIPForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2);
            AddComponent(elasticTension);

            var pressure = new PressureGradientForm(SpeciesName, d);
            AddComponent(pressure);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => material.Density;

        public override string CodomainName => codomainName;

    }
}

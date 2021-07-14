using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class LinearNavierCauchy : BulkEquation {

        string speciesName;

        Solid material;

        string codomainName;

        public LinearNavierCauchy(string speciesName, Solid material, int d, int D, double artificialViscosity) {
            this.speciesName = speciesName;
            this.material = material;
            this.codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);

            var convection = new LinearConvectionForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], D, material.Density);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D)[d]);
            AddParameter(BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D)[d]);
            AddComponent(convection);

            if(d == 0) {
                var elasticTension = new LinearIncompressibleNeoHookeanX(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), material.Lame2);
                AddComponent(elasticTension);
            } else if(d == 1) {
                var elasticTension = new LinearIncompressibleNeoHookeanY(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), material.Lame2);
                AddComponent(elasticTension);
            } else {
                throw new Exception("Spatial Dimension not supported.");
            }

            var pressure = new PressureGradientForm(SpeciesName, d);
            AddComponent(pressure);

            var artificialViscosityForm = new SIPForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, artificialViscosity);
            AddComponent(artificialViscosityForm);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => material.Density;

        public override string CodomainName => codomainName;

    }
}


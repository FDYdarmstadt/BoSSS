using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class NavierCauchy : BulkEquation {

        string speciesName;

        Solid material;

        string codomainName;

        public NavierCauchy(string speciesName, Solid material, int d, int D) {
            this.speciesName = speciesName;
            this.material = material;
            this.codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            AddVariableNames(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.Pressure);
            
            var convection = new NonLinearConvectionForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), D, material.Density);
            AddComponent(convection);

            //var elasticTension = new IncompressibleNeoHookean(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D),d, material.Lame2);
            //*
            if (d == 0) {
                var elasticTension = new LinearIncompressibleNeoHookeanX(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), material.Lame2);
                AddComponent(elasticTension);
            } else if (d == 1) {
                var elasticTension = new LinearIncompressibleNeoHookeanY(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), material.Lame2);
                AddComponent(elasticTension);
            } else {
                throw new Exception("Spatial Dimension not supported.");
            }
            //*/

            var pressure = new PressureGradientForm(SpeciesName, d);
            AddComponent(pressure);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => material.Density;

        public override string CodomainName => codomainName;

    }
}

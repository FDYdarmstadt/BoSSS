using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    /// <summary>
    /// Momentum equation for the Solid /full linearization to be used with <see cref="BoSSS.Solution.Control.NonLinearSolverCode.Newton"/>
    /// </summary>
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

            
            var convection = new NonLinearConvectionForm(SpeciesName, 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], 
                BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D),
                d, material.Density);
            AddComponent(convection);
            
            var pressure = new PressureGradientForm(SpeciesName, d);
            AddComponent(pressure);

            var eulerAlmansi0 = new SIPForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2);
            AddComponent(eulerAlmansi0);

            var eulerAlmansi1 = new SIPTransposeForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2);
            AddComponent(eulerAlmansi1);
            
            //if(d == 0) {
            var viscosity = new SIPForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, material.Viscosity);
            AddComponent(viscosity);
            //} else {
            //    var viscosity = new Fake_ipFlux(1.3, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], material.Viscosity);
            //    AddComponent(viscosity);
            //}
            string gravity = BoSSS.Solution.NSECommon.VariableNames.GravityVector(D)[d];
            string gravityOfSpecies = gravity + "#" + SpeciesName;
            var gravityComponent = new BoSSS.Solution.XNSECommon.Operator.MultiPhaseSource(gravityOfSpecies, speciesName);
            AddComponent(gravityComponent);
            AddParameter(gravityOfSpecies);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => material.Density;

        public override string CodomainName => codomainName;

    }
}

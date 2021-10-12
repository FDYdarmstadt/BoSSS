using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    /// <summary>
    /// Momentum equation for the Solid / ad-hoc linearization to be used with <see cref="BoSSS.Solution.Control.NonLinearSolverCode.Picard"/>
    /// </summary>
    class LinearNavierCauchy : BulkEquation {

        string speciesName;

        Solid material;

        string codomainName;

        
        public LinearNavierCauchy(string speciesName, Solid material, int d, int D) {
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
            //Console.WriteLine("##################### Rem: nix convection.");
            
            

            var pressure = new PressureGradientForm(SpeciesName, d);
            AddComponent(pressure);
            /*
            // laplacian of displacement:
            var eulerAlmansi0 = new SIPForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2);
            //eulerAlmansi0.PenaltySafety = 0.0;
            AddComponent(eulerAlmansi0);

            var eulerAlmansi1 = new SIPTransposeForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2);
            AddComponent(eulerAlmansi1);
            //eulerAlmansi1.PenaltySafety = 0.0;
            */Console.WriteLine("##################### Rem: displacement coupling deakt.");

            var viscosity = new SIPForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, material.Viscosity);
            AddComponent(viscosity);

            //var velocityBoundaryPenalty = new EdgePenaltyForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[d], 1);
            //AddComponent(velocityBoundaryPenalty);

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


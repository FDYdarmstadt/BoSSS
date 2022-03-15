using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.XNSECommon;
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

        static public double EulerAlamansiPenalty = 1.0;

        public NavierCauchy(string speciesName, Solid material, int d, int D, IncompressibleMultiphaseBoundaryCondMap boundaryMap) {
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
            
            var pressure = new GradientForm(SpeciesName, d, BoSSS.Solution.NSECommon.VariableNames.Pressure);
            AddComponent(pressure);

            if(material.Lame2 != 0.0)
            {
                var eulerAlmansi0 = new SIPForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2, EulerAlamansiPenalty);
                AddComponent(eulerAlmansi0);

                var eulerAlmansi1 = new SIPTransposeForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2, EulerAlamansiPenalty);
                AddComponent(eulerAlmansi1);
                
                /*
                var gradUGradUT = new GradATGradBForm(SpeciesName, ZwoLevelSetSolver.VariableNames.DisplacementVector(D),
                    ZwoLevelSetSolver.VariableNames.DisplacementVector(D), d, material.Lame2, EulerAlamansiPenalty);
                AddComponent(gradUGradUT);
                */
                
            }
            if(material.Viscosity != 0.0)
            {
                var viscosity = new SIPForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, material.Viscosity);
                AddComponent(viscosity);
                var viscosityT = new SIPTransposeForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, material.Viscosity);
                AddComponent(viscosityT);
            }
            
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

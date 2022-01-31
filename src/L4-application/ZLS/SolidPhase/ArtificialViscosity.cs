using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    
    
    class ArtificialViscosity : BulkEquation {

        string speciesName;

        //Solid material;

        string codomainName;

        public ArtificialViscosity(string speciesName, double viscosity, int d, int D) {
            this.speciesName = speciesName;
            //this.material = material;
            this.codomainName = BoSSS.Solution.NSECommon.EquationNames.MomentumEquationComponent(d);
            AddVariableNames(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            
            var artificialViscosity = new SIPForm(SpeciesName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D), d, viscosity);
            AddComponent(artificialViscosity);
        }

        public override string SpeciesName => speciesName;

        public override double MassScale => 1.0;// material.Density;

        public override string CodomainName => codomainName;

    }
}

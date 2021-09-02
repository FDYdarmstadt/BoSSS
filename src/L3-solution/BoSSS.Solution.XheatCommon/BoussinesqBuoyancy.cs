using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XheatCommon {
    
    /// <summary>
    /// Buoyancy approximation in otherwise incompressible flow
    /// </summary>
    public class BoussinesqApproximation_Buoyancy : IVolumeForm, ISupportsJacobianComponent, ISpeciesFilter {
        string species;
        string[] ParameterName;
        public BoussinesqApproximation_Buoyancy(string species, string parameter, int d, ThermalParameters thermparams) {
            this.species = species;
            this.d = d;
            ParameterName = new string[] { parameter };
            switch (species) {
                case "A": rho = thermparams.rho_A; alpha = thermparams.alpha_A; T0 = thermparams.T_sat; break;
                case "B": rho = thermparams.rho_B; alpha = thermparams.alpha_B; T0 = thermparams.T_sat; ; break;
                default: throw new ArgumentException("Unknown species.");
            }
        }
        double rho, alpha, T0;
        int d;

        public string ValidSpecies => species;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new string[] { BoSSS.Solution.NSECommon.VariableNames.Temperature };

        public IList<string> ParameterOrdering => ParameterName;

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return -alpha / rho * cpv.Parameters[0] * (U[0] - T0) * V;
        }
    }
    
}

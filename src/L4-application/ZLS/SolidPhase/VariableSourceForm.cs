using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    /// <summary>
    /// Minimal volume source term
    /// </summary>
    public class VariableSourceForm : IVolumeForm, ISupportsJacobianComponent {
        string[] variableName;
        double scale;

        public VariableSourceForm(string variable, double scale = 1) {
            variableName = new string[]
            {
                variable
            };
            this.scale = scale;
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering => variableName;

        public IList<string> ParameterOrdering => new string[0];

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return scale * U[0] * V;
        }
    }

    public class MultiPhaseVariableSource : VariableSourceForm, ISpeciesFilter {
        string species;

        public MultiPhaseVariableSource(string species, string variable, double scale = 1) : base(variable, scale) {
            this.species = species;
        }

        public string ValidSpecies => species;
    }
}

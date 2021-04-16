using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XNSECommon.Operator {


    /// <summary>
    /// Minimal volume source term
    /// </summary>
    public class ParameterSource : IVolumeForm, ISupportsJacobianComponent {
        string[] parameterName;
        double scale;

        public ParameterSource(string parameter, double scale = 1.0) {
            parameterName = new string[]
            {
                parameter
            };
            this.scale = scale;
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ArgumentOrdering => new string[0];

        public IList<string> ParameterOrdering => parameterName;

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return cpv.Parameters[0] * V * scale;
        }
    }

    public class MultiPhaseSource : ParameterSource, ISpeciesFilter {
        string species;

        public MultiPhaseSource(string parameter, string species, double scale = 1.0) : base(parameter, scale) {
            this.species = species;
        }

        public string ValidSpecies => species;
    }
}

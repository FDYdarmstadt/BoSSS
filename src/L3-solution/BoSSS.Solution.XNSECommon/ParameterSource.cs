using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XNSECommon.Operator
{
    public class ParameterSource : IVolumeForm
    {
        string[] parameterName;

        public ParameterSource(string parameter)
        {
            parameterName = new string[]
            {
                parameter
            };
        }

        public TermActivationFlags VolTerms
        {
            get
            {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ArgumentOrdering => new string[0];

        public IList<string> ParameterOrdering => parameterName;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            return cpv.Parameters[0] * V;
        }
    }

    public class MultiPhaseSource : ParameterSource, ISpeciesFilter
    {
        string species;

        public MultiPhaseSource(string parameter, string species) : base(parameter)
        {
            this.species = species;
        }

        public string ValidSpecies => species;
    }
}

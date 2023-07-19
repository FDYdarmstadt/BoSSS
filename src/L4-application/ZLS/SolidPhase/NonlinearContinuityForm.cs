using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;

namespace ZwoLevelSetSolver.SolidPhase {
    class NonlinearContinuityForm : IVolumeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;
        string[] variables;
        int d;

        public NonlinearContinuityForm(string speciesName, string[] variableNames) {
            this.speciesName = speciesName;
            this.variables = variableNames;
            this.d = d;
        }

        public IList<string> ArgumentOrdering => variables;

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV; }
        }

        public string ValidSpecies => speciesName;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double flux = GradU[0, 0] * GradU[1,1] - GradU[0, 1] * GradU[1, 0];
            return flux * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] {DivergenceDerivVol };
        }
    }
}

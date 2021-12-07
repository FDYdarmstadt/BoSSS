using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class DivergenceGradientForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent {
        string species;
        int d;
        int D;

        public DivergenceGradientForm(string species, int d, int D) {
            this.species = species;
            this.d = d;
            this.D = D;
        }

        public IList<string> ArgumentOrdering => VariableNames.DisplacementVector(D);

        public IList<string> ParameterOrdering => new string[0];

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.GradUxV; }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0 * inp.Normal[d] * _vA;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double divU = 0;
            for(int i = 0; i < D; ++i) {
                divU += 0.5 * _Grad_uIN[i, i];
                divU += 0.5 * _Grad_uOUT[i, i];
            }
            return divU * inp.Normal[d] * (_vIN - _vOUT);
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxGradV; }
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double divU = 0;
            for(int i = 0; i < D; ++i) {
                divU += GradU[i, i];
            }
            double acc = divU * GradV[d];
            return -acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}

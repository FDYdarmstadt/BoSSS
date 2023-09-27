using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class PressureGradientForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent {
        string species;
        int d;

        public PressureGradientForm(string species, int d) {
            this.species = species;
            this.d = d;
        }

        public TermActivationFlags VolTerms { 
            get { return TermActivationFlags.UxGradV; } 
        }

        public IList<string> ArgumentOrdering => new string[] { BoSSS.Solution.NSECommon.VariableNames.Pressure };

        public IList<string> ParameterOrdering => new string[0];

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV|TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV|TermActivationFlags.V; }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return _uA[0] * inp.Normal[d] * _vA;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return (_uIN[0] + _uOUT[0]) * 0.5 * inp.Normal[d] * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = U[0] * GradV[d];
            return -acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}

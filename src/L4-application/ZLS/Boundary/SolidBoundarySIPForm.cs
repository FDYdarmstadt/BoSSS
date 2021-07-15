using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.Boundary {

    class SolidBoundarySIPForm : ILevelSetForm {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        int D;
        double viscosity;

        public SolidBoundarySIPForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double viscosity) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.D = D;
            
            variableNames = VariableNames.DisplacementVector(D);
            this.viscosity = viscosity;
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.GradUxV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc1 = 0.0;
            for(int i = 0; i < D; i++) {
                //acc1 -= 0.5 * viscosity * (_Grad_uOUT[d, i]) * ( _vIN - _vOUT) * inp.Normal[i];  // consistency term  
                acc1 -= viscosity * (_Grad_uOUT[d, i]) * ( _vOUT) * inp.Normal[i];  // consistency term  
                //acc1 -= 0.5 * viscosity * (_Grad_vOUT[i]) * (_uIN[d] - _uOUT[d]) * inp.Normal[i];  // symmetry term
            }
            return acc1;
        }
    }


}

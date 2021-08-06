using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class SIPForm : IVolumeForm, IEdgeForm, ISpeciesFilter, IEquationComponentCoefficient {
        double viscosity;
        string species;
        int d;
        string[] variableNames;

        public SIPForm(string species, string[] variables, int d, double viscosity) {
            this.species = species;
            this.viscosity = viscosity;
            this.variableNames = variables;
            this.d = d;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.AllOn; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.GradUxV| TermActivationFlags.UxGradV| TermActivationFlags.V |TermActivationFlags.GradV);}
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV); }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double acc1 = 0.0;
            double pnlty = 2 * Penalty(inp.jCellIn, -1);

            Vector dirichlet = new Vector(D);

            for (int i = 0; i < D; i++) {
                acc1 -= viscosity * _Grad_uIN[d, i] * _vIN  * inp.Normal[i];  // consistency term  
                acc1 -= viscosity * _Grad_vIN[i] * (_uIN[d] - dirichlet[d]) * inp.Normal[i];  // symmetry term
            }
            acc1 += (_uIN[d] - dirichlet[d]) * _vIN  * pnlty * viscosity;
            return acc1;
        }

        MultidimensionalArray cj;

        double penalty;

        int D;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = cs.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice
            
            cj = cs.CellLengthScales;
        }

        double Penalty(int jCellIn, int jCellOut) {
            double penaltySizeFactor_A = 1.0 / cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc1 = 0.0;
            for(int i = 0; i < D; i++) {
                acc1 -= 0.5 * viscosity * ( _Grad_uIN[d, i] + _Grad_uOUT[d, i]) * (_vIN - _vOUT) * inp.Normal[i];  // consistency term  
                acc1 -= 0.5 * viscosity * (_Grad_vIN[i] + _Grad_vOUT[i]) * (_uIN[d] - _uOUT[d]) * inp.Normal[i];  // symmetry term
            }

            double pnlty = Penalty(inp.jCellIn, inp.jCellOut);
            acc1 += (_uIN[d] - _uOUT[d]) * (_vIN - _vOUT) * pnlty * viscosity;
            return acc1;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int i = 0; i < D; ++i) {
                acc += GradU[d, i] * GradV[i];
            }
            acc *= viscosity;
            return acc ;
        }
    }
}

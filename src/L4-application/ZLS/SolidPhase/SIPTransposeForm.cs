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
    public class SIPTransposeForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {
        double viscosity;
        string species;
        int d;
        string[] variableNames;
        public double PenaltySafety;
        public SIPTransposeForm(string species, string[] variables, int d, double viscosity, double __PenaltySafety = 4.0) {
            this.PenaltySafety = __PenaltySafety;
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
            get { return TermActivationFlags.AllOn; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.AllOn; }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double acc1 = 0.0;
            for (int i = 0; i < D; i++) {
                acc1 -=  viscosity * (_Grad_uIN[i, d] ) * (_vIN ) * inp.Normal[i];  // consistency term  
                acc1 -=  viscosity * (_Grad_vIN[i] ) * (_uIN[i] - 0) * inp.Normal[d];  // symmetry term
            }

            double pnlty = PenaltyIn(inp.jCellIn);
            acc1 += PenaltySafety * (_uIN[d] - 0) * (_vIN) * pnlty * viscosity;
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

        double PenaltyIn(int jCellIn) {
            double penaltySizeFactor_A = 1.0 / cj[jCellIn];

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_A * penalty;
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        double PenaltyOut(int jCellOut) {
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_B * penalty;
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc1 = 0.0;
            for (int i = 0; i < D; i++) {
                acc1 -= 0.5 * viscosity * (_Grad_uIN[i, d] + _Grad_uOUT[i, d]) * (_vIN - _vOUT) * inp.Normal[i];  // consistency term  
                acc1 -= 0.5 * viscosity * (_Grad_vIN[i] + _Grad_vOUT[i]) * (_uIN[i] - _uOUT[i]) * inp.Normal[d];  // symmetry term
            }
            double penalty = Math.Max(PenaltyIn(inp.jCellIn), PenaltyOut(inp.jCellOut));
            acc1 += PenaltySafety * penalty * (_uIN[d] - _uOUT[d]) * (_vIN - _vOUT) * viscosity;
            return acc1;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int i = 0; i < D; ++i) {
                acc += GradU[i, d] * GradV[i];
            }
            acc *= viscosity;
            return acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}

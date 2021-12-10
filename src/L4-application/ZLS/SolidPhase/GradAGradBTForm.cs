using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace ZwoLevelSetSolver.SolidPhase {
    class GradAGradBTForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {
        protected double viscosity;
        string species;
        protected int d;
        protected string[] variableNames;
        public double PenaltySafety;

        public GradAGradBTForm(string species, string[] variablesA, string[] variablesB, int d, double viscosity, double __PenaltySafety = 4.0) {
            this.species = species;
            this.viscosity = viscosity;
            this.PenaltySafety = __PenaltySafety;
            this.variableNames = variablesA.Cat(variablesB);
            this.d = d;
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.V | TermActivationFlags.GradV); }
        }

        public string ValidSpecies => species;

        public virtual double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double acc1 = 0.0;
            for(int i = 0; i < D; i++) {
                double GradUGradU = 0;
                for(int j = 0; j < D; ++j) {
                    GradUGradU += 0.5 * _Grad_uIN[d, j] * _Grad_uIN[D + i, j];
                }
                acc1 -= viscosity * GradUGradU * (_vIN) * inp.Normal[i];  // consistency term  
            }
            //double penalty = Math.Max(PenaltyIn(inp.jCellIn), -1);
            //acc1 += (_uIN[d] + _uIN[D+d]) * PenaltySafety * penalty * (_vIN) * viscosity;
            return acc1;
        }

        MultidimensionalArray cj;

        double penalty;

        protected int D;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = cs.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            this.cj = cs.CellLengthScales;
            //this.cj = ((BoSSS.Foundation.Grid.Classic.GridData)(cs.GrdDat)).Cells.CellLengthScale;
        }

        protected double PenaltyIn(int jCellIn) {
            double penaltySizeFactor_A = 1.0 / cj[jCellIn];

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_A * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        double PenaltyOut(int jCellOut) {
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_B * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV); }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc1 = 0.0;
            for(int i = 0; i < D; i++) {
                double GradUGradU = 0;
                for(int j = 0; j < D; ++j) {
                    GradUGradU += 0.5 * _Grad_uIN[d, j] * _Grad_uIN[D + i, j];
                    GradUGradU += 0.5 * _Grad_uOUT[d, j] * _Grad_uOUT[D + i, j];
                }
                acc1 -= viscosity * GradUGradU * (_vIN - _vOUT) * inp.Normal[i];  // consistency term  
            }
            //double penalty = Math.Max(PenaltyIn(inp.jCellIn), PenaltyOut(inp.jCellOut));
            //acc1 += (_uIN[d] - _uOUT[d]) * PenaltySafety * penalty * (_vIN - _vOUT) * viscosity;
            //acc1 += (_uIN[D+d] - _uOUT[D+d]) * PenaltySafety * penalty * (_vIN - _vOUT) * viscosity;
            return acc1;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxGradV; }
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int i = 0; i < D; ++i) {
                double GradUGradU = 0;
                for(int j = 0; j < D; ++j) {
                    GradUGradU += GradU[d, j] * GradU[D + i, j];
                }
                acc += GradUGradU * GradV[i];
            }
            acc *= viscosity;
            return acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }
    }
}

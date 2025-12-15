using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HFSISolver.SolidPhase {

    class EdgePenaltyForm : IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {

        string speciesName;
        string[] variables;
        double scale = 1;

        public EdgePenaltyForm(string speciesName, string variableName, double scale = 1) {
            this.speciesName = speciesName;
            this.variables = new string[] { variableName };
            if(scale.IsNaNorInf()) {
                throw new ArgumentException("Illegal scaling value for EdgePenaltyForm: " + scale);
            }
            this.scale = scale;
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV; }
        }

        public IList<string> ArgumentOrdering => variables;

        public IList<string> ParameterOrdering => null;

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
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
            double penaltySizeFactor = cj[jCellIn];
            if(jCellOut > -1) {
                penaltySizeFactor = Math.Max(penaltySizeFactor, cj[jCellOut]);
            }
            Debug.Assert(!double.IsNaN(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor * penalty;
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double flux =  (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT);
            flux *= scale * Penalty(inp.jCellIn, inp.jCellOut);
            Debug.Assert(!flux.IsNaNorInf());
            return flux;

        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class EdgePenaltyForm1 : IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {

        string speciesName;
        string[] variables;
        double scale = 1;

        public EdgePenaltyForm1(string speciesName, string variableName, double scale = 1) {
            this.speciesName = speciesName;
            this.variables = new string[] { variableName };
            this.scale = scale;
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV; }
        }

        public IList<string> ArgumentOrdering => variables;

        public IList<string> ParameterOrdering => null;

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
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
            double penaltySizeFactor = 1/ cj[jCellIn];
            if(jCellOut > -1) {
                penaltySizeFactor = Math.Max(penaltySizeFactor, 1/cj[jCellOut]);
            }
            Debug.Assert(!double.IsNaN(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double flux = (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT);
            flux *= scale * Penalty(inp.jCellIn, inp.jCellOut);
            return flux;

        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}

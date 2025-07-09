using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    class EdgePenaltyForm : IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {

        string speciesName;
        string[] variables;
        double scale = 1;

        public EdgePenaltyForm(string speciesName, string variableName, double scale = 1) {
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
            //return scale * Penalty(inp.jCellIn, -1) * (_uA[0]) * (_vA);
            //return scale * (_uA[0]) * (_vA);

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
            double penaltySizeFactor_A = 1/cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 && !cj[jCellOut].IsNaNorInf() ? 1/cj[jCellOut] : 0;
            // 0.0 case:
            // there is no OUT-cell, i.e., current edge is at the boundary of the domain
            // -- or -- 
            // in 3D, with certain cut-cell quad rules (e.g. Algoim),
            // when the Level Set passes exactly through the corner of the cell,
            // it might happen that the species-volume of the OUT-cell is empty,
            // but the edge integral is only almost zero (slightly positive, e.g., weights around 10e-14 or so).
            // In such cases, the edge integral is evaluated (since it is a non-empty rule),
            // but the OUT-cell already has a NAN-length scale assigned.

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor;
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double flux =  (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT);
            flux *= scale * Penalty(inp.jCellIn, inp.jCellOut);
            return flux;

        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}

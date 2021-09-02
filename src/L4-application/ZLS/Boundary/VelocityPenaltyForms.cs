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
    class NoSlipVelocityPenaltyForm : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        double viscosity;
        double lame2;

        public NoSlipVelocityPenaltyForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.viscosity = viscosity;
            this.lame2 = lame2;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.UxV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        MultidimensionalArray cj;

        double penalty;

        int D;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = csB.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = csB.CellLengthScales;
        }

        double Penalty(int jCellIn, int jCellOut) {
            double penaltySizeFactor_A = 1/cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1/cj[jCellOut] : 0;

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return 200;
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double flux =  (_uIN[d] - _uOUT[d]) * (_vIN - _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            flux *= Math.Max(viscosity, lame2);
            //double flux = (viscosity * (_uIN[d] - _uOUT[d]) * _vIN - lame2 * (_uIN[d] - _uOUT[d]) * _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            //double flux =  (_uIN[d] - _uOUT[d]) * (_vIN ) * Penalty(inp.jCellIn, inp.jCellOut);
            return flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class SlipVelocityPenaltyForm : ILevelSetForm, ILevelSetEquationComponentCoefficient {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        double viscosity;
        double lame2;

        public SlipVelocityPenaltyForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.viscosity = viscosity;
            this.lame2 = lame2;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.UxV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        MultidimensionalArray cj;

        double penalty;

        int D;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = csB.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = csB.CellLengthScales;
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
            double flux = 0;
            for(int i = 0; i < D; ++i) {
                flux += (_uIN[i] - _uOUT[i]) * inp.Normal[i];
            }
            flux *= inp.Normal[d] * (_vIN - _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            flux *= Math.Max(viscosity, lame2);
            //double flux = (viscosity * (_uIN[d] - _uOUT[d]) * _vIN - lame2 * (_uIN[d] - _uOUT[d]) * _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            //double flux =  (_uIN[d] - _uOUT[d]) * (_vIN ) * Penalty(inp.jCellIn, inp.jCellOut);
            return flux;
        }
    }

    class NavierSlipVelocityPenaltyForm : ILevelSetForm, ILevelSetEquationComponentCoefficient {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        double viscosity;
        double lame2;
        double slipLength;

        public NavierSlipVelocityPenaltyForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double viscosity, double lame2, double slipLength) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.viscosity = viscosity;
            this.lame2 = lame2;
            this.slipLength = slipLength;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.UxV|TermActivationFlags.GradUxV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        MultidimensionalArray cj;

        double penalty;

        int D;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = csB.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = csB.CellLengthScales;
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
            Vector normalVelocityGradient = new Vector(D);
            for(int i = 0; i < D; ++i) {
                for(int j = 0; j < D; ++j) {
                    normalVelocityGradient[j] += (_Grad_uIN[j, i]) * (-inp.Normal[i]);
                }
            }
            double slip = normalVelocityGradient[d];
            for(int i = 0; i < D; ++i) {
                slip -= (inp.Normal[d]) * (inp.Normal[i]) * normalVelocityGradient[i];
            }
            slip *= slipLength;
            double flux = (_uIN[d] - slip - _uOUT[d]) * (_vIN - _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            flux *= Math.Max(viscosity, lame2);
            //double flux = (viscosity * (_uIN[d] - _uOUT[d]) * _vIN - lame2 * (_uIN[d] - _uOUT[d]) * _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            //double flux =  (_uIN[d] - _uOUT[d]) * (_vIN ) * Penalty(inp.jCellIn, inp.jCellOut);
            return flux;
        }
    }
}

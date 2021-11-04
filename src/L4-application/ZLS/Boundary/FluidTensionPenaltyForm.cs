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
    class FluidTensionPenaltyForm : ILevelSetForm, ILevelSetEquationComponentCoefficient {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        double viscosity;
        double lame2;

        public FluidTensionPenaltyForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.viscosity = viscosity;
            this.lame2 = lame2;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
            variableNames = variableNames.Cat(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
            variableNames = variableNames.Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure);
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.GradUxV; }
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
            double tension = 0;

            // ((fluidTension - solidTension) * n 
            for(int i = 0; i < D; i++) {
                // = consistency terms
                tension -= lame2 * (_Grad_uOUT[D + d, i]) * inp.Normal[i];
                tension += viscosity * (_Grad_uIN[d, i]) * inp.Normal[i];
            }
            //pressure penalty
            tension += 1 * (_uOUT[2 * D]) * inp.Normal[d];
            tension -= 1 * (_uIN[2 * D]) * inp.Normal[d];

            //Test functions and penalty

            tension *= 0.001 * Penalty(inp.jCellIn, inp.jCellOut) * (-_vIN);
            return tension;
        }
    }
}


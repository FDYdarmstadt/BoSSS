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

namespace HFSISolver.Boundary {
    class NeoHookeanBoundaryForm : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D;//Should check here!!!
        int d;
        double fluidViscosity;
        double solidViscosity;
        double lame2;

        public NeoHookeanBoundaryForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, 
            double fluidViscosity, double solidViscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.solidViscosity = solidViscosity;
            this.lame2 = lame2;
            this.fluidViscosity = fluidViscosity;
            this.d = d;
            this.D = D;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
            variableNames = variableNames.Cat(HFSISolver.VariableNames.DisplacementVector(D));
            variableNames = variableNames.Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure);
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            //pressure
            double pIn = (_uIN[2 * D]) * inp.Normal[d];
            double pOut = (_uOUT[2 * D]) * inp.Normal[d];

            //Tension, consistency
            double fluidStress = 0.0;
            double fluidStressT = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 1 * fluidViscosity * (_Grad_uIN[d, i]) * inp.Normal[i];
                fluidStressT -= 1 * fluidViscosity * (_Grad_uIN[i, d]) * inp.Normal[i];
            }

            double solidStress = 0.0;
            double solidStressT = 0.0;
            double viscousStress = 0.0;
            double viscousStressT = 0.0;

            for(int i = 0; i < D; i++) {
                solidStress -= 1 * lame2 * (_Grad_uOUT[D + d, i]) * inp.Normal[i];
                solidStressT -= 1 * lame2 * (_Grad_uOUT[D + i, d]) * inp.Normal[i];
                
                viscousStress -= 1 * solidViscosity * (_Grad_uOUT[d, i]) * inp.Normal[i];
                viscousStressT -= 1 * solidViscosity * (_Grad_uOUT[i, d]) * inp.Normal[i];
            }

            double stress = (fluidStress + fluidStressT + pIn + solidStress + solidStressT + viscousStress + viscousStressT + pOut) * 0.5;//Lacking surface tension and curvature here?

            return  (stress ) * _vIN -  (stress) * _vOUT;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class NeoHookeanNeumannForm : ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        int d;
        double fluidViscosity;
        double solidViscosity;
        double lame2;

        public NeoHookeanNeumannForm(string fluidSpecies, string solidSpecies, int d, int levelSetIndex, 
            double fluidViscosity, double solidViscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.solidViscosity = solidViscosity;
            this.lame2 = lame2;
            this.fluidViscosity = fluidViscosity;
            this.d = d;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
            variableNames = variableNames.Cat(HFSISolver.VariableNames.DisplacementVector(D));
            variableNames = variableNames.Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure);
        }

        public int LevelSetIndex => levelSetIndex;

        public string PositiveSpecies => solidSpecies;

        public string NegativeSpecies => fluidSpecies;

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV; }
        }

        MultidimensionalArray cj;

        double penalty;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
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

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            //pressure
            double pIn = (_uIN[2 * D]) * inp.Normal[d];

            //Tension, consistency
            double fluidStress = 0.0;
            double fluidStressT = 0.0;

            for(int i = 0; i < D; i++) {
                fluidStress -= 1 * fluidViscosity * (_Grad_uIN[d, i]) * inp.Normal[i];
                fluidStressT -= 1 * fluidViscosity * (_Grad_uIN[i, d]) * inp.Normal[i];
            }

            return (fluidStress + pIn ) * (_vIN)
                + ( fluidStress + pIn) * (-_vOUT);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

}

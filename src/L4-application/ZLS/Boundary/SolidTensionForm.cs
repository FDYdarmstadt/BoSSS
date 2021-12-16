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
using NSEVariableNames = BoSSS.Solution.NSECommon.VariableNames;
namespace ZwoLevelSetSolver.Boundary {
    class SolidTensionForm : ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        double artificialViscosity;
        double fluidViscosity;
        double solidViscosity;

        public SolidTensionForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double artificialViscosity, double fluidViscosity, double solidViscosity) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.solidViscosity = solidViscosity;
            this.fluidViscosity = fluidViscosity;

            this.d = d;
            this.artificialViscosity = artificialViscosity;
            this.variableNames = VariableNames.DisplacementVector(D).Cat(NSEVariableNames.VelocityVector(D)).Cat(NSEVariableNames.Pressure);
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
            //pressure
            double pIn = (_uIN[2 * D]) * inp.Normal[d];
            double pOut = (_uOUT[2 * D]) * inp.Normal[d];

            //Tension, consistency
            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 1 * fluidViscosity * (_Grad_uIN[D+d, i]) * inp.Normal[i];
            }
            double viscousStress = 0.0;

            for(int i = 0; i < D; i++) {
                viscousStress -= 1 * solidViscosity * (_Grad_uOUT[D+d, i]) * inp.Normal[i];
                viscousStress -= 1 * solidViscosity * (_Grad_uOUT[D+i, d]) * inp.Normal[i];
            }

            // ((fluidTension + solidTension) * n 
            for(int i = 0; i < D; i++) {
                // = consistency terms
                //tension -= viscosity * _Grad_uOUT[d, i] * inp.Normal[i];
                tension -= artificialViscosity * _Grad_uOUT[d, i] * inp.Normal[i];
            }
            //tension += viscosity * Penalty(inp.jCellIn, inp.jCellOut) * (_uIN[d] - _uOUT[d]);

            //Impose equality by weakly imposing via boundary conditions: solidStress = fluidStress + pIn - viscousStress - pOut
            //tension = artificialViscosity * (-pIn - fluidStress + pOut + viscousStress);

            return tension * (-_vOUT);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class ExtensionSolidTensionForm : ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int d;
        double extensionViscosity;
        double artificialViscosity;

        public ExtensionSolidTensionForm(string fluidSpecies, string solidSpecies, string[] variableNames, int d, int D, int levelSetIndex, double extensionViscosity, double artificialViscosity) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.extensionViscosity = extensionViscosity;
            this.artificialViscosity = artificialViscosity;
            this.variableNames = variableNames;
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
            double extensionTension = 0;
            double artificialTension = 0;
            // ((fluidTension + solidTension) * n 
            for(int i = 0; i < D; i++) {
                // = consistency terms
                extensionTension -= extensionViscosity * _Grad_uIN[d, i] * inp.Normal[i];
                artificialTension -= artificialViscosity * _Grad_uOUT[d, i] * inp.Normal[i];
            }
            double penalty = 0;
            penalty += Math.Max(extensionViscosity, artificialViscosity) * Penalty(inp.jCellIn, inp.jCellOut) * (_uIN[d] - _uOUT[d]) * (_vIN); ; 

            return extensionTension * _vIN - artificialTension * _vOUT + penalty;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}


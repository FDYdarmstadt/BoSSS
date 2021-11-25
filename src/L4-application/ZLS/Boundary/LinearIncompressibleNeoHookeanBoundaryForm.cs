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
    class SolidLinearIncompressibleNeoHookeanBoundaryForm : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        int d;
        double fluidViscosity;
        double solidViscosity;
        double lame2;

        public SolidLinearIncompressibleNeoHookeanBoundaryForm(string fluidSpecies, string solidSpecies, int d, int levelSetIndex, double fluidViscosity, double solidViscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.solidViscosity = solidViscosity;
            this.lame2 = lame2;
            this.fluidViscosity = fluidViscosity;
            this.d = d;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
            variableNames = variableNames.Cat(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
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
            double pressure = 0.0;

            //pressure
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[d];
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[d];

            //Tension, consistency
            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 0.5 * fluidViscosity * (_Grad_uIN[d, i]) * inp.Normal[i];
            }

            double solidStress = 0.0;
            for (int i = 0; i < D; i++) {
                solidStress -= 0.5 * lame2 * (_Grad_uOUT[D + d, i]) * inp.Normal[i];
                solidStress -= 0.5 * lame2 * (_Grad_uOUT[D + i, d]) * inp.Normal[i];

                solidStress -= 0.5 * solidViscosity * (_Grad_uOUT[d, i]) * inp.Normal[i];
            }

            return (pressure + fluidStress + solidStress) * (_vIN - _vOUT);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class SlipSolidLinearIncompressibleNeoHookeanBoundaryForm : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        int d;
        double fluidViscosity;
        double solidViscosity;
        double lame2;

        public SlipSolidLinearIncompressibleNeoHookeanBoundaryForm(string fluidSpecies, string solidSpecies, int d, int levelSetIndex, double fluidViscosity, double solidViscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.solidViscosity = solidViscosity;
            this.lame2 = lame2;
            this.fluidViscosity = fluidViscosity;
            this.d = d;
            variableNames = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D);
            variableNames = variableNames.Cat(ZwoLevelSetSolver.VariableNames.DisplacementVector(D));
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
            double pressure = 0.0;

            //pressure
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[d];
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[d];

            //Tension, consistency
            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                for(int j = 0; j < D; ++j) {
                    fluidStress -= 0.5 * fluidViscosity *inp.Normal[i] * (_Grad_uIN[i, j]) * inp.Normal[j] * inp.Normal[d];
                }
            }

            double solidStress = 0.0;
            for(int i = 0; i < D; i++) {
                for(int j = 0; j < D; ++j) {
                    solidStress -= 0.5 * lame2 * inp.Normal[i] * (_Grad_uOUT[D + i, j]) * inp.Normal[j] * inp.Normal[d];
                    solidStress -= 0.5 * lame2 * inp.Normal[i] * (_Grad_uOUT[D + j, i]) * inp.Normal[j] * inp.Normal[d];

                    solidStress -= 0.5 * solidViscosity * inp.Normal[i] * (_Grad_uOUT[i, j]) * inp.Normal[j] * inp.Normal[d];
                }
            }

            return (pressure + fluidStress + solidStress) * (_vIN - _vOUT);
        }


        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}

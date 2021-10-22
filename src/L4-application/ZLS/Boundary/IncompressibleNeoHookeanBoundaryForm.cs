using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.Boundary {
    class IncompressibleNeoHookeanBoundaryForm : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        int d;
        double viscosity;
        double lame2;

        LeftCauchyGreenDeformationTensor deformation;
        IncompressibleNeoHookeanSolid cauchyStress;

        MultidimensionalArray outerStress;
        MultidimensionalArray outerDeformation;


        public IncompressibleNeoHookeanBoundaryForm(string fluidSpecies, string solidSpecies, int levelSetIndex, int d, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.lame2 = lame2;
            this.viscosity = viscosity;
            variableNames = ZwoLevelSetSolver.VariableNames.DisplacementVector(D);
            variableNames = variableNames.Cat(BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
            variableNames = variableNames.Cat(BoSSS.Solution.NSECommon.VariableNames.Pressure);

            cauchyStress = new IncompressibleNeoHookeanSolid(lame2, D);
            deformation = new LeftCauchyGreenDeformationTensor(D);

            outerStress = MultidimensionalArray.Create(D, D);
            outerDeformation = MultidimensionalArray.Create(D, D);
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
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[d];
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[d];

            //Tension, consistency
            double fluidStress = 0.0;
            for (int i = 0; i < D; i++) {
                fluidStress -= 0.5 * viscosity * (_Grad_uIN[D + d, i]) * inp.Normal[i];
            }
            
            double solidStress = 0.0;
            deformation.Calculate(_Grad_uOUT, outerDeformation);
            cauchyStress.Calculate(outerDeformation, outerStress);

            for (int i = 0; i < D; i++) {
                solidStress = -0.5 * ( outerStress[d, i]) * inp.Normal[i];
            }

            return (pressure + fluidStress + solidStress) * (_vIN -_vOUT);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }
}

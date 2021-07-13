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
    class SolidLinearIncompressibleNeoHookeanBoundaryFormX : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        double viscosity;
        double lame2;

        public SolidLinearIncompressibleNeoHookeanBoundaryFormX(string fluidSpecies, string solidSpecies,int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.lame2 = lame2;
            this.viscosity = viscosity;
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
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[0];
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[0];

            //Tension, consistency
            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 0.5 * viscosity * (_Grad_uIN[0, i]) * inp.Normal[i];
            }
            double solidStress = 0.0;
            solidStress += 0.5 * lame2 * (- _Grad_uOUT[D + 1, 1] * inp.Normal[0] +  _Grad_uOUT[D + 0, 1] * inp.Normal[1]);
            solidStress += 0.5 * lame2 * (-_Grad_uOUT[D + 1, 1] * inp.Normal[0] + _Grad_uOUT[D + 1, 0] * inp.Normal[1]);

            solidStress *= -1;

            //Tension, symmetry
            
            /*
            for(int i = 0; i < D; i++) {
                acc2 -= 0.5 * viscosity * (_Grad_vOUT[i]) * inp.Normal[i];
            }
            acc2 *= -_uIN[0];
            //*/
            /*
            acc2 -= 1 * lame2 * ( _Grad_vOUT[1]) * inp.Normal[1] *  _uOUT[D + 0];
            acc2 += 1 * lame2 * ( _Grad_vOUT[0]) * inp.Normal[1] *  _uOUT[D + 1];
            acc2 -= 1 * lame2 * ( _Grad_vOUT[1] * inp.Normal[0] - _Grad_vOUT[0] * inp.Normal[1]) * _uOUT[D + 1];
            //*/
            return (pressure + fluidStress+ solidStress) * -_vOUT;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class SolidLinearIncompressibleNeoHookeanBoundaryFormY : ILevelSetForm, ISupportsJacobianComponent{
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        double viscosity;
        double lame2;

        public SolidLinearIncompressibleNeoHookeanBoundaryFormY(string fluidSpecies, string solidSpecies, int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.lame2 = lame2;
            this.viscosity = viscosity;
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
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[1];
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[1];

            //Tension, consistency
            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 0.5 * viscosity * (_Grad_uIN[1, i]) * inp.Normal[i];
            }
            double solidStress = 0.0;
            solidStress += 0.5 * lame2 * (( _Grad_uOUT[D + 1, 0]) * inp.Normal[0] - ( _Grad_uOUT[D + 0, 0]) * inp.Normal[1]);
            solidStress += 0.5 * lame2 * (_Grad_uOUT[D + 0, 1] * inp.Normal[0] - _Grad_uOUT[D + 0, 0] * inp.Normal[1]);
            
            solidStress *= -1;
            
            //Tension, symmetry
            
            /*
            for(int i = 0; i < D; i++) {
                acc2 -= 0.5 * viscosity * (_Grad_vOUT[i]) * inp.Normal[i];
            }
            acc2 *= -_uIN[1];
            //*/
            /*
            acc2 += 1 * lame2 * (_Grad_vOUT[1]) * inp.Normal[0] * _uOUT[D + 0];
            acc2 -= 1 * lame2 * (_Grad_vOUT[0]) * inp.Normal[0] * _uOUT[D + 1];
            acc2 -= 1 * lame2 * (-_Grad_vOUT[1] * inp.Normal[0] +  _Grad_vOUT[0] * inp.Normal[1]) * _uOUT[D + 0];
            //*/
            return (pressure + fluidStress+ solidStress) * -_vOUT;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class FluidLinearIncompressibleNeoHookeanBoundaryFormX : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        double viscosity;
        double lame2;

        public FluidLinearIncompressibleNeoHookeanBoundaryFormX(string fluidSpecies, string solidSpecies, int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
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
            get { return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];

        

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double pressure = 0;
            //pressure
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[0];
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[0];


            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 0.5 * viscosity * (_Grad_uIN[0, i]) * inp.Normal[i];
            }

            // Tension, consistency
            double solidStress = 0.0;
            
            solidStress += 0.5 * lame2 * (-_Grad_uOUT[D + 1, 1] * inp.Normal[0] + _Grad_uOUT[D + 0, 1] * inp.Normal[1]);
            solidStress += 0.5 * lame2 * (-_Grad_uOUT[D + 1, 1] * inp.Normal[0] + _Grad_uOUT[D + 1, 0] * inp.Normal[1]);

            solidStress *= -1;

            //Tension, symmetry
            double fluidSymmetry = 0.0;
            /*
            for(int i = 0; i < D; i++) {
                fluidSymmetry -= 1 * viscosity * (_Grad_vIN[i]) * inp.Normal[i];
            }
            fluidSymmetry *= _uIN[0] - _uOUT[0];
            
            /*
            acc1 -= 0.5 * (_Grad_vIN[1]) * inp.Normal[1] * (_uOUT[D + 0]);
            acc1 += 0.5 * (_Grad_vIN[0]) * inp.Normal[1] * (_uOUT[D + 1]);
            acc1 -= 0.5 * ((_Grad_vIN[1] ) * inp.Normal[0] - (_Grad_vIN[0]) * inp.Normal[1]) * (_uOUT[D + 1]);

            //*/
            return (pressure + fluidStress + solidStress) * _vIN + fluidSymmetry;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class FluidLinearIncompressibleNeoHookeanBoundaryFormY : ILevelSetForm, ISupportsJacobianComponent {
        int levelSetIndex;
        string solidSpecies;
        string fluidSpecies;
        string[] variableNames;
        int D = 2;
        double viscosity;
        double lame2;

        public FluidLinearIncompressibleNeoHookeanBoundaryFormY(string fluidSpecies, string solidSpecies, int levelSetIndex, double viscosity, double lame2) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
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
            get { return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[0];



        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double pressure = 0.0;

            //pressure
            pressure += 0.5 * (_uOUT[2 * D]) * inp.Normal[1];
            pressure += 0.5 * (_uIN[2 * D]) * inp.Normal[1];

            //Tension, consistency
            double fluidStress = 0.0;
            for(int i = 0; i < D; i++) {
                fluidStress -= 0.5 * viscosity * (_Grad_uIN[1, i]) * inp.Normal[i];
            }

            double solidStress = 0.0;
            solidStress += 0.5 * lame2 * ((_Grad_uOUT[D + 1, 0]) * inp.Normal[0] - (_Grad_uOUT[D + 0, 0]) * inp.Normal[1]);
            solidStress += 0.5 * lame2 * (_Grad_uOUT[D + 0, 1] * inp.Normal[0] - _Grad_uOUT[D + 0, 0] * inp.Normal[1]);

            solidStress *= -1;

            //Tension, symmetry
            double fluidSymmetry = 0.0;
            /*
            for(int i = 0; i < D; i++) {
                fluidSymmetry -= 1* viscosity * (_Grad_vIN[i]) * inp.Normal[i];
            }
            fluidSymmetry *=  _uIN[1]- _uOUT[1];
            
            /*
            acc2 += 0.5 * (_Grad_vIN[1]) * inp.Normal[0] * (_uOUT[D + 0]);
            acc2 -= 0.5 * (_Grad_vIN[0]) * inp.Normal[0] * (_uOUT[D+ 1]);
            acc2 -= 0.5 * (-(_Grad_vIN[1]) * inp.Normal[0] + (_Grad_vIN[0]) * inp.Normal[1]) * (_uOUT[D + 0]);

            //*/
            return (pressure + fluidStress + solidStress) *_vIN + fluidSymmetry;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}

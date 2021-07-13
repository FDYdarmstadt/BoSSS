using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class LinearIncompressibleNeoHookeanX : IVolumeForm, IEdgeForm, ISpeciesFilter, IEquationComponentCoefficient, ISupportsJacobianComponent {
        double viscosity;
        string species;
        int D = 2;
        string[] variableNames;

        public LinearIncompressibleNeoHookeanX(string species, string[] variables, double viscosity) {
            this.species = species;
            this.viscosity = viscosity;
            this.variableNames = variables;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV); }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV); }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double acc1 = 0.0;
            acc1 -= (-(_Grad_uIN[1, 1]) * inp.Normal[0] + (_Grad_uIN[1, 0]) * inp.Normal[1]) * _vIN;
            acc1 -= (_Grad_vIN[1]) * inp.Normal[1] * (_uIN[0]);
            acc1 += (_Grad_vIN[0]) * inp.Normal[1] * (_uIN[1]);

            acc1 -=  (-(_Grad_uIN[1, 1] ) * inp.Normal[0] + (_Grad_uIN[1, 0]) * inp.Normal[1]) * (_vIN);
            acc1 -=  ((_Grad_vIN[1]) * inp.Normal[0] - (_Grad_vIN[0]) * inp.Normal[1]) * (_uIN[1]);

            acc1 *= -1;
            double pnlty = Penalty(inp.jCellIn, -1);
            acc1 += (_uIN[0]) * (_vIN) * pnlty;
            return viscosity * acc1;

        }

        MultidimensionalArray cj;

        double penalty;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            Debug.Assert(cs.GrdDat.SpatialDimension == 2, "Only 2d support");
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;
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
            double acc1 = 0.0;

            //adj(grad U)
            acc1 += 0.5 * (-(_Grad_uIN[1, 1] + _Grad_uOUT[1,1]) * inp.Normal[0] + (_Grad_uIN[0, 1] + _Grad_uOUT[0,1])* inp.Normal[1] ) * (_vIN - _vOUT);
            acc1 += 0.5 * (_Grad_vIN[1] + _Grad_vOUT[1]) * inp.Normal[1] * (_uIN[0] -_uOUT[0]); 
            acc1 += -0.5 * (_Grad_vIN[0] + _Grad_vOUT[0]) * inp.Normal[1] * (_uIN[1] - _uOUT[1]);

            //adj(grad U)^T
            acc1 += 0.5 * (-(_Grad_uIN[1, 1] + _Grad_uOUT[1, 1]) * inp.Normal[0] + (_Grad_uIN[1, 0] + _Grad_uOUT[1, 0]) * inp.Normal[1]) * (_vIN - _vOUT);
            acc1 += 0.5 *( (_Grad_vIN[1] + _Grad_vOUT[1]) * inp.Normal[0] - (_Grad_vIN[0] + _Grad_vOUT[0]) * inp.Normal[1]) * (_uIN[1] - _uOUT[1]);

            acc1 *= -1;

            double pnlty = 4 * Penalty(inp.jCellIn, inp.jCellOut);
            acc1 += (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT) * pnlty;
            return viscosity* acc1;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            //adj(grad U)
            double acc = -GradU[1, 1] * GradV[0] + GradU[0, 1] * GradV[1];

            //adj(grad U)^T
            acc += -GradU[1, 1] * GradV[0] + GradU[1, 0] * GradV[1];

            acc *= viscosity;
            return acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class LinearIncompressibleNeoHookeanY : IVolumeForm, IEdgeForm, ISpeciesFilter, IEquationComponentCoefficient, ISupportsJacobianComponent {
        double viscosity;
        string species;
        int D = 2;
        string[] variableNames;

        public LinearIncompressibleNeoHookeanY(string species, string[] variables, double viscosity) {
            this.species = species;
            this.viscosity = viscosity;
            this.variableNames = variables;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV); }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV); }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double acc1 = 0.0;
            acc1 -= ((_Grad_uIN[1, 0] ) * inp.Normal[0] - (_Grad_uIN[0, 0]) * inp.Normal[1]) * (_vIN);
            acc1 += (_Grad_vIN[1] ) * inp.Normal[0] * (_uIN[0]);
            acc1 -= (_Grad_vIN[0] ) * inp.Normal[0] * (_uIN[1]);

            acc1 -= ((_Grad_uIN[0, 1] ) * inp.Normal[0] - (_Grad_uIN[0, 0] ) * inp.Normal[1]) * (_vIN );
            acc1 -= (-(_Grad_vIN[1] ) * inp.Normal[0] + (_Grad_vIN[0] ) * inp.Normal[1]) * (_uIN[0]);
            acc1 *= -1;

            double pnlty = Penalty(inp.jCellIn, -1);
            acc1 += (_uIN[1]) * (_vIN) *  pnlty ;
            return acc1 * viscosity;

        }

        MultidimensionalArray cj;

        double penalty;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            Debug.Assert(cs.GrdDat.SpatialDimension == 2, "Only 2d support");
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;
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
            double acc1 = 0.0;

            //adj(grad U)
            acc1 += 0.5 * ((_Grad_uIN[1, 0] + _Grad_uOUT[1, 0]) * inp.Normal[0] - (_Grad_uIN[0, 0] + _Grad_uOUT[0, 0]) * inp.Normal[1]) * (_vIN - _vOUT);
            acc1 += -0.5 * (_Grad_vIN[1] + _Grad_vOUT[1]) * inp.Normal[0] * (_uIN[0] - _uOUT[0]);
            acc1 += 0.5 * (_Grad_vIN[0] + _Grad_vOUT[0]) * inp.Normal[0] * (_uIN[1] - _uOUT[1]);

            //adj(grad U)^T
            acc1 += 0.5 * ((_Grad_uIN[0, 1] + _Grad_uOUT[0, 1]) * inp.Normal[0] - (_Grad_uIN[0, 0] + _Grad_uOUT[0, 0]) * inp.Normal[1]) * (_vIN - _vOUT);
            acc1 += 0.5 * (-(_Grad_vIN[1] + _Grad_vOUT[1]) * inp.Normal[0] + (_Grad_vIN[0] + _Grad_vOUT[0]) * inp.Normal[1]) * (_uIN[0] - _uOUT[0]);

            acc1 *= -1;

            double pnlty = 4 * Penalty(inp.jCellIn, inp.jCellOut);
            acc1 += (_uIN[1] - _uOUT[1]) * (_vIN - _vOUT) * pnlty;
            return acc1 * viscosity;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            //adj(grad U)
            double acc = GradU[1, 0] * GradV[0] - GradU[0, 0] * GradV[1];
            
            //adj(grad U)^T
            acc += GradU[0, 1] * GradV[0] - GradU[0, 0] * GradV[1];

            acc *= viscosity;
            return acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}


using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
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
        double penaltySolid;
        double penaltyFluid;

        public NoSlipVelocityPenaltyForm(string fluidSpecies, string solidSpecies, int d, int D, int levelSetIndex, double penaltyFluid, double penaltySolid) {
            this.levelSetIndex = levelSetIndex;
            this.fluidSpecies = fluidSpecies;
            this.solidSpecies = solidSpecies;
            this.d = d;
            this.penaltyFluid = penaltyFluid;
            this.penaltySolid = penaltySolid;
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

        MultidimensionalArray cjOut;

        MultidimensionalArray cjIn;

        double penalty;

        int D;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = csB.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cjOut = csB.CellLengthScales;
            cjIn = csA.CellLengthScales;
        }

        double PenaltyIn(int jCellIn) {
            double penaltySizeFactor_A = 1/cjIn[jCellIn];
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_A * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        double PenaltyOut(int jCellOut) {
            double penaltySizeFactor_B = jCellOut >= 0 ? 1 / cjOut[jCellOut] : 0;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_B * penalty;
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            //return 200;
            return µ;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, 
            double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            double flux =  PenaltyIn(inp.jCellIn) * (_uIN[d] - _uOUT[d])  * (penaltyFluid * _vIN - penaltySolid *_vOUT);
            return flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class SlipVelocityPenaltyForm : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {
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

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

    class NavierSlipVelocityPenaltyForm : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {
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
            //Vector normalVelocityGradient = new Vector(D);
            //for(int i = 0; i < D; ++i) {
            //    for(int j = 0; j < D; ++j) {
            //        normalVelocityGradient[j] += (_Grad_uIN[j, i]) * (-inp.Normal[i]);
            //        //normalVelocityGradient[j] += (_Grad_uIN[j, i] + _Grad_uIN[i, j]) * (-inp.Normal[i]);
            //    }
            //}
            //double slip = normalVelocityGradient[d];
            //for(int i = 0; i < D; ++i) {
            //    slip -= (inp.Normal[d]) * (inp.Normal[i]) * normalVelocityGradient[i];
            //}
            double[] grad_u_n_tangential = ComputeGradUNormalTangentialProjection(_Grad_uIN, inp.Normal);
            double slip = grad_u_n_tangential[d];
            slip *= slipLength;

            //double uIN_tangential = _uIN[d];
            //double uOUT_tangential = _uOUT[d];
            //for(int i = 0; i < D; ++i) {
            //    uIN_tangential -= (uIN_tangential * inp.Normal[i]) * inp.Normal[i];
            //    uOUT_tangential -= (uIN_tangential * inp.Normal[i]) * inp.Normal[i];
            //}
            double[] uIN_tangential = ComputeTangentialVelocity(_uIN, inp.Normal);
            double[] uOUT_tangential = ComputeTangentialVelocity(_uOUT, inp.Normal);

            double flux = (uIN_tangential[d] + slip - uOUT_tangential[d]) * (_vIN - _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            //double flux = (_uIN[d] - slip - _uOUT[d]) * (_vIN - _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            flux *= viscosity;
            //flux *= Math.Max(viscosity, lame2);
            //double flux = (viscosity * (_uIN[d] - _uOUT[d]) * _vIN - lame2 * (_uIN[d] - _uOUT[d]) * _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);
            //double flux =  (_uIN[d] - _uOUT[d]) * (_vIN ) * Penalty(inp.jCellIn, inp.jCellOut);

            //// Enforce normal velocity continuity
            //double uIN_normal = 0.0, uOUT_normal = 0.0;
            //for(int i = 0; i < D; ++i) {
            //    uIN_normal += _uIN[d] * inp.Normal[i] * inp.Normal[i];
            //    uOUT_normal += _uOUT[d] * inp.Normal[i] * inp.Normal[i];
            //}
            double[] uIN_normal = ComputeNormalVelocity(_uIN, inp.Normal);
            double[] uOUT_normal = ComputeNormalVelocity(_uOUT, inp.Normal);

            flux += viscosity * (uIN_normal[d] - uOUT_normal[d]) * (_vIN - _vOUT) * Penalty(inp.jCellIn, inp.jCellOut);

            return flux;
        }

        double[] ComputeTangentialVelocity(double[] _uIN, double[] normal) {
            int D = _uIN.Length;
            double[] u_tangential = new double[D];

            //  I - n ⊗ n
            double[,] P = new double[D, D];

            for(int i = 0; i < D; i++) {
                for(int j = 0; j < D; j++) {
                    if(i == j)
                        P[i, j] = 1.0 - normal[i] * normal[j];
                    else
                        P[i, j] = -normal[i] * normal[j];
                }
            }

            // u_tangential = P * _uIN
            for(int i = 0; i < D; i++) {
                u_tangential[i] = 0.0;
                for(int j = 0; j < D; j++) {
                    u_tangential[i] += P[i, j] * _uIN[j];
                }
            }

            return u_tangential;
        }

        double[] ComputeNormalVelocity(double[] _uIN, double[] normal) {
            int D = _uIN.Length;
            double[] u_normal = new double[D];

            // Compute dot product (u · n)
            double dot = 0.0;
            for(int i = 0; i < D; ++i) {
                dot += _uIN[i] * normal[i];
            }

            // u_normal = (u · n) * n
            for(int i = 0; i < D; ++i) {
                u_normal[i] = dot * normal[i];
            }

            return u_normal;
        }

        double[] ComputeGradUNormalTangentialProjection(double[,] grad_u, double[] normal) {
            int D = normal.Length;
            double[] grad_u_n = new double[D];         // ∇u ⋅ n
            double[] tangential = new double[D];       // final tangential projection

            // Step 1: Compute grad_u_n = grad_u ⋅ n
            for(int i = 0; i < D; ++i) {
                grad_u_n[i] = 0.0;
                for(int j = 0; j < D; ++j) {
                    grad_u_n[i] += (grad_u[i, j] + grad_u[j, i]) * normal[j];
                }
            }

            // Step 2: Compute tangential = (I - n⊗n) ⋅ grad_u_n
            // Equivalent to: tangential[i] = grad_u_n[i] - (n ⋅ grad_u_n) * n[i]

            double dot = 0.0;
            for(int i = 0; i < D; ++i)
                dot += grad_u_n[i] * normal[i];

            for(int i = 0; i < D; ++i)
                tangential[i] = grad_u_n[i] - dot * normal[i];

            return tangential;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var JacobiComp = new LevelSetFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { JacobiComp };
        }
    }

}

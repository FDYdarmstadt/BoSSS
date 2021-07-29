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
    class IncompressibleNeoHookean : IVolumeForm, IEdgeForm, ISpeciesFilter, IEquationComponentCoefficient, ISupportsJacobianComponent {

        double lame2;
        string species;
        int D = 2;
        int d;
        string[] variableNames;

        LeftCauchyGreenDeformationTensor deformation;
        IncompressibleNeoHookeanSolid cauchyStress;

        MultidimensionalArray innerStress;
        MultidimensionalArray outerStress;
        MultidimensionalArray innerDeformation;
        MultidimensionalArray outerDeformation;

        public IncompressibleNeoHookean(string species, string[] variables, int d, double lame2) {
            this.species = species;
            this.lame2 = lame2;
            this.variableNames = variables;
            this.d = d;

            cauchyStress = new IncompressibleNeoHookeanSolid(lame2, D);
            deformation = new LeftCauchyGreenDeformationTensor(D);

            innerStress = MultidimensionalArray.Create(D, D);
            outerStress = MultidimensionalArray.Create(D, D);
            innerDeformation = MultidimensionalArray.Create(D, D);
            outerDeformation = MultidimensionalArray.Create(D, D);
            gradVIn = new double[D, D];
            gradVOut = new double[D, D];
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

            //Consistency
            deformation.Calculate(_Grad_uIN, innerDeformation);
            cauchyStress.Calculate(innerDeformation, innerStress);

            for (int i = 0; i < D; i++) {
                acc1 = - (innerStress[d, i]) * (_vIN) * inp.Normal[i];
            }

            //Symmetry
            Embedd(_Grad_vIN, gradVIn, d);
            deformation.Calculate(gradVIn, innerDeformation);
            cauchyStress.Calculate(innerDeformation, innerStress);

            for (int i = 0; i < D; i++) {
                acc1 = - (innerStress[d, i]) * (_uIN[d]) * inp.Normal[i];
            }
            double pnlty = Penalty(inp.jCellIn, -1);
            acc1 += _uIN[d] * _vIN * pnlty * lame2;
            return acc1;
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
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        double[,] gradVIn;
        double[,] gradVOut;

        static void Embedd(double[] gradV, double[,] result, int d) {
            for(int i = 0; i < gradV.Length; ++i) {
                result[d, i] = gradV[i];
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc1 = 0.0;

            //Consistency
            deformation.Calculate(_Grad_uIN, innerDeformation);
            deformation.Calculate(_Grad_uOUT, outerDeformation);

            cauchyStress.Calculate(innerDeformation, innerStress);
            cauchyStress.Calculate(outerDeformation, outerStress);

            for (int i = 0; i < D; ++i) {
                acc1 = -0.5 * (innerStress[d,i] + outerStress[d,i]) * (_vIN - _vOUT) * inp.Normal[i];
            }

            /*
            //Symmetry
            Embedd(_Grad_vIN, gradVIn, d);
            Embedd(_Grad_vOUT, gradVOut, d);

            
            //deformation.Calculate(gradVIn, innerDeformation);
            //deformation.Calculate(gradVOut, outerDeformation);

            cauchyStress.Calculate(innerDeformation, innerStress);
            cauchyStress.Calculate(outerDeformation, outerStress);

            for (int i = 0; i < D; ++i) {
                acc1 = -0.5 * (innerStress[d, i] + outerStress[d, i]) * (_uIN[d] - _uOUT[d]) * inp.Normal[i];
            }
            */

            double pnlty = Penalty(inp.jCellIn, inp.jCellOut);
            acc1 += lame2 * (_uIN[d] - _uOUT[d]) * (_vIN - _vOUT) * pnlty;

            Debug.Assert(!double.IsNaN(acc1));
            return acc1;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            deformation.Calculate(GradU, innerDeformation);
            cauchyStress.Calculate(innerDeformation, innerStress);

            double acc = 0;
            for (int i = 0; i < D; ++i) {
                acc += innerStress[d, i] * GradV[i];
            }

            return acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }
    }

    class LeftCauchyGreenDeformationTensor {

        int D;

        MultidimensionalArray buffer;


        public LeftCauchyGreenDeformationTensor(int D) {
            this.D = D;
            buffer = MultidimensionalArray.Create(D, D);
        }

        public void Calculate(double[,] gradDisplacement, MultidimensionalArray result) {
            // Tensor= (1 - gradU)^-1
            DeformationGradient(gradDisplacement, D, result);
            result.InvertInPlace();
            
            // buffer= (1 - gradU)^-T
            buffer.CopyFrom(result);
            buffer.TransposeInPlace();

            //Tensor = (1 - gradU)^-1 * (1 - gradU)^-T
            result = result *  buffer;
        }

        static void DeformationGradient(double[,] gradDisplacement, int D, MultidimensionalArray result) {
            result.Clear();
            for(int i = 0; i < D; ++i) {
                result[i, i] = 1;
                for (int j = 0; j < D; ++j) {
                    result[i, j] -= gradDisplacement[i, j];
                }
            }
        }
    }

    /// <summary>
    /// Without pressure
    /// </summary>
    class IncompressibleNeoHookeanSolid {

        double lame2;

        public IncompressibleNeoHookeanSolid(double lame2, int D) {
            this.lame2 = lame2;
        }

        public void Calculate(MultidimensionalArray leftCauchyGreenDeformationTensor, MultidimensionalArray result) {
            result.CopyFrom(leftCauchyGreenDeformationTensor);
            result.Scale(2 * lame2);
        }
    }
}

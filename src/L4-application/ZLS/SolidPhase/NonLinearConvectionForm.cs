using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    //All boundaries are a Wall for now
    class NonLinearConvectionForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;

        double rho;

        string[] variableNames;

        string[] parameternames;

        int D;

        public NonLinearConvectionForm(string speciesName, string variableName, string[] velocity, int D, double rho) {
            this.speciesName = speciesName;
            this.variableNames = new string[] { variableName }.Cat(velocity);
            this.D = D;
            this.rho = rho;
            this.parameternames = new string[] { };
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.UxGradV | TermActivationFlags.GradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => parameternames;

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {

            double r = 0.0;
            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += _uIN[0] * (_uIN[1+0] * inp.Normal[0] + _uIN[1 + 1] * inp.Normal[1]);
            if(D == 3) {
                r += _uIN[0] * _uIN[1 + 2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            Vector VelocityMeanIn = new Vector(D);
            for(int d = 0; d < D; ++d) {
                VelocityMeanIn[d] = _uIN[1 + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double uJump = _uIN[0];
            r += LambdaIn * uJump;

            r *= rho;
            return r * (_vIN);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double r = 0.0;

            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += _uIN[0] * (_uIN[1+0] * inp.Normal[0] + _uIN[1 + 1] * inp.Normal[1]);
            r += _uOUT[0] * (_uOUT[1+0] * inp.Normal[0] + _uOUT[1+1] * inp.Normal[1]);
            if(D == 3) {
                r += _uIN[0] * _uIN[1+2] * inp.Normal[2] + _uOUT[0] * _uOUT[1+2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            Vector VelocityMeanIn = new Vector(D);
            Vector VelocityMeanOut = new Vector(D);
            for(int d = 0; d < D; ++d) {
                VelocityMeanIn[d] = _uIN[1 + d];
                VelocityMeanOut[d] = _uOUT[1 + d];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal);
            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = _uIN[0] - _uOUT[0];
            r += Lambda * uJump;
            r *= 0.5;

            r *= rho;
            return r * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < D; d++)
                acc += U[0] * U[1+d] * GradV[d];
            acc *= rho;
            return -acc;
        }
    }
}

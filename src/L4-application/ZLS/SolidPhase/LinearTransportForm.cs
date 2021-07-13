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
    public class LinearTransportForm : IVolumeForm, IEdgeForm, ISpeciesFilter {

        string speciesName;

        double rho;

        string[] variableNames;

        string[] parameternames;

        int D;

        int d;

        public LinearTransportForm(string speciesName, string[] variableNames, int d, int D, double rho) {
            this.speciesName = speciesName;
            this.variableNames = variableNames;
            this.D = D;
            this.d = d;
            this.rho = rho;
            this.parameternames = BoSSS.Solution.NSECommon.VariableNames.Velocity0Vector(D).Cat(
                BoSSS.Solution.NSECommon.VariableNames.Velocity0MeanVector(D));
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
            r += inp.Parameters_IN[d] * (_uIN[0] * inp.Normal[0] + _uIN[1] * inp.Normal[1]);
            if(D == 3) {
                r += inp.Parameters_IN[d] * _uIN[2] * inp.Normal[2];
            }
            // Calculate dissipative part
            // ==========================

            Vector VelocityMeanIn = new Vector(D);
            for(int i = 0; i < D; ++i) {
                VelocityMeanIn[i] = inp.Parameters_IN[D + i];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double uJump = _uIN[d];
            r += LambdaIn * uJump;

            r *= rho;
            return r * (_vIN);
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double r = 0.0;

            r += inp.Parameters_IN[d] * (_uIN[0] * inp.Normal[0] + _uIN[1] * inp.Normal[1]);
            r += inp.Parameters_OUT[d] * (_uOUT[0] * inp.Normal[0] + _uOUT[1] * inp.Normal[1]);
            if(D == 3) {
                r += inp.Parameters_IN[d] * _uIN[2] * inp.Normal[2] + inp.Parameters_OUT[d] * _uOUT[2] * inp.Normal[2];
            }

            // Calculate dissipative part
            // ==========================

            Vector VelocityMeanIn = new Vector(D);
            Vector VelocityMeanOut = new Vector(D);
            for(int i = 0; i < D; ++i) {
                VelocityMeanIn[i] = inp.Parameters_IN[D + i];
                VelocityMeanOut[i] = inp.Parameters_OUT[D + i];
            }

            double LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normal);
            double LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normal);
            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double uJump = _uIN[d] - _uOUT[d];
            //r += Lambda * uJump;

            r *= rho;
            r *= 0.5;
            return r * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int i = 0; i < D; i++)
                acc += U[i] * GradV[i];
            acc *= cpv.Parameters[d];
            acc *= rho;
            return -acc;
        }
    }
}
